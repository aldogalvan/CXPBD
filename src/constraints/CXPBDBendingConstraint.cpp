//
// Created by aldof on 4/30/2022.
//

#include "CXPBDBendingConstraint.h"
#include <algorithm>

inline Eigen::Matrix3d convertVecToCrossOp(const Eigen::Vector3d& vec)
{
    Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();

    mat(0, 1) = -vec(2);
    mat(0, 2) = +vec(1);
    mat(1, 0) = +vec(2);
    mat(1, 2) = -vec(0);
    mat(2, 0) = -vec(1);
    mat(2, 1) = +vec(0);

    return mat;
};

inline double calculateCotTheta(const Eigen::Vector3d& x, const Eigen::Vector3d& y)
{
    const double scaled_cos_theta = x.dot(y);
    const double scaled_sin_theta = x.cross(y).norm();
    return scaled_cos_theta / scaled_sin_theta;
}

cXPBDBendingConstraint::scalar_type
cXPBDBendingConstraint::evaluate(positions_type const& p, masses_type const& M) const
{
    const Eigen::Vector3d& x_0 = p.row(indices()[0]);
    const Eigen::Vector3d& x_1 = p.row(indices()[1]);
    const Eigen::Vector3d& x_2 = p.row(indices()[2]);
    const Eigen::Vector3d& x_3 = p.row(indices()[3]);

    const Eigen::Vector3d p_10 = x_1 - x_0;
    const Eigen::Vector3d p_20 = x_2 - x_0;
    const Eigen::Vector3d p_30 = x_3 - x_0;

    const Eigen::Vector3d n_0 = p_10.cross(p_20).normalized();
    const Eigen::Vector3d n_1 = p_10.cross(p_30).normalized();

    double clamp = n_0.dot(n_1);
    if(clamp < -1.0)
    {
        clamp = -1.0;
    }
    else if (clamp > 1.0)
    {
        clamp = 1.0;
    }

    const double current_dihedral_angle = std::acos(clamp);

    assert(n_0.norm() > 0.0);
    assert(n_1.norm() > 0.0);
    assert(!std::isnan(current_dihedral_angle));

    return current_dihedral_angle - m_di_ang;
}

void cXPBDBendingConstraint::project(
        positions_type& p,
        positions_type& p0_,
        masses_type const& M,
        scalar_type& lagrange,
        scalar_type const dt,
        gradient_type& F) const
{
    const Eigen::Vector3d& x_0 = p.row(indices()[0]);
    const Eigen::Vector3d& x_1 = p.row(indices()[1]);
    const Eigen::Vector3d& x_2 = p.row(indices()[2]);
    const Eigen::Vector3d& x_3 = p.row(indices()[3]);

    // Assuming that p_0 = [ 0, 0, 0 ]^T without loss of generality
    const Eigen::Vector3d p_1 = x_1 - x_0;
    const Eigen::Vector3d p_2 = x_2 - x_0;
    const Eigen::Vector3d p_3 = x_3 - x_0;

    const Eigen::Vector3d p_1_cross_p_2 = p_1.cross(p_2);
    const Eigen::Vector3d p_1_cross_p_3 = p_1.cross(p_3);

    const Eigen::Vector3d n_0 = p_1_cross_p_2.normalized();
    const Eigen::Vector3d n_1 = p_1_cross_p_3.normalized();

    const double d = n_0.dot(n_1);

    Eigen::MatrixXd grad_C(4,3);

    // If the current dihedral angle is sufficiently small or large (i.e., zero or pi), return zeros.
    // This is only an ad-hoc solution for stability and it needs to be solved in a more theoretically grounded way.
    constexpr double epsilon = 1e-12;
    if (1.0 - d * d < epsilon)
    {
        grad_C.setZero();
        return;
    }

    const double common_coeff = -1.0 / std::sqrt(1.0 - d * d);

    auto calc_grad_of_normalized_cross_prod_wrt_p_a =
            [](const Eigen::Vector3d& p_a, const Eigen::Vector3d& p_b, const Eigen::Vector3d& n) -> Eigen::Matrix3d
            {
                return +(1.0 / p_a.cross(p_b).norm()) * (-convertVecToCrossOp(p_b) + n * (n.cross(p_b)).transpose());
            };

    auto calc_grad_of_normalized_cross_prod_wrt_p_b =
            [](const Eigen::Vector3d& p_a, const Eigen::Vector3d& p_b, const Eigen::Vector3d& n) -> Eigen::Matrix3d
            {
                return -(1.0 / p_a.cross(p_b).norm()) * (-convertVecToCrossOp(p_a) + n * (n.cross(p_a)).transpose());
            };

    const Eigen::Matrix3d partial_n_0_per_partial_p_1 = calc_grad_of_normalized_cross_prod_wrt_p_a(p_1, p_2, n_0);
    const Eigen::Matrix3d partial_n_1_per_partial_p_1 = calc_grad_of_normalized_cross_prod_wrt_p_a(p_1, p_3, n_1);
    const Eigen::Matrix3d partial_n_0_per_partial_p_2 = calc_grad_of_normalized_cross_prod_wrt_p_b(p_1, p_2, n_0);
    const Eigen::Matrix3d partial_n_1_per_partial_p_3 = calc_grad_of_normalized_cross_prod_wrt_p_b(p_1, p_3, n_1);

    const Eigen::Vector3d grad_C_wrt_p_1 =
            common_coeff * (partial_n_0_per_partial_p_1.transpose() * n_1 + partial_n_1_per_partial_p_1.transpose() * n_0);
    const Eigen::Vector3d grad_C_wrt_p_2 = common_coeff * partial_n_0_per_partial_p_2.transpose() * n_1;
    const Eigen::Vector3d grad_C_wrt_p_3 = common_coeff * partial_n_1_per_partial_p_3.transpose() * n_0;
    const Eigen::Vector3d grad_C_wrt_p_0 = -grad_C_wrt_p_1 - grad_C_wrt_p_2 - grad_C_wrt_p_3;

    grad_C.row(0) = grad_C_wrt_p_0;
    grad_C.row(1) = grad_C_wrt_p_1;
    grad_C.row(2) = grad_C_wrt_p_2;
    grad_C.row(3) = -grad_C_wrt_p_3;

    // Calculate the constraint function value
    const double C = evaluate(p,M);


    // Skip if the gradient is sufficiently small
    constexpr double very_small_value = 1e-12;
    if (grad_C.norm() < very_small_value)
    {
        return;
    }

    // Calculate time-scaled compliance
    const double alpha_tilde = alpha_ / (dt * dt);

    double w0 = M(indices()[0]);
    double w1 = M(indices()[1]);
    double w2 = M(indices()[2]);
    double w3 = M(indices()[3]);

    //
    double weighted_sum_of_gradients = w0*grad_C_wrt_p_0.squaredNorm() + w1*grad_C_wrt_p_1.squaredNorm() +
            w2*grad_C_wrt_p_2.squaredNorm() + w3*grad_C_wrt_p_3.squaredNorm();

    // Calculate \Delta lagrange multiplier
    const double delta_lagrange =
            (-C - alpha_tilde * lagrange) /
            (weighted_sum_of_gradients + alpha_tilde);

    // Calculate \Delta x
    p.row(indices()[0]) += w0*delta_lagrange*grad_C_wrt_p_0;
    p.row(indices()[1]) += w1*delta_lagrange*grad_C_wrt_p_1;
    p.row(indices()[2]) += w2*delta_lagrange*grad_C_wrt_p_2;
    p.row(indices()[3]) += w3*delta_lagrange*grad_C_wrt_p_3;


    // Update the lagrange multiplier
    lagrange += delta_lagrange;

}