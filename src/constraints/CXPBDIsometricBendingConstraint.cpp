//
// Created by aldo on 6/16/22.
//

#include "CXPBDIsometricBendingConstraint.h"

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

cXPBDIsometricBendingConstraint::scalar_type
cXPBDIsometricBendingConstraint::evaluate(positions_type const& p, masses_type const& M) const
{
    double sum = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
        for (unsigned int j = 0; j < 4; ++j)
        {
            sum += m_Q(i, j) * double(p.row(indices()[i]) * p.row(indices()[j]).transpose());
        }
    }
    return 0.5 * sum;
}

void
cXPBDIsometricBendingConstraint::project(positions_type& V, positions_type& V0, masses_type const& M,
        scalar_type& lagrange, scalar_type const dt, gradient_type& F) const
{


    Eigen::MatrixXd grad_C(4,3);

    for (unsigned int i = 0; i < 4; ++i)
    {
        Eigen::Vector3d sum = Eigen::Vector3d::Zero();
        for (unsigned int j = 0; j < 4; ++j)
        {
            sum += m_Q(i, j) * V.row(indices()[j]);
        }
        grad_C.row(i) = sum;
    }

    // Skip if the gradient is sufficiently small
    constexpr double very_small_value = 1e-12;
    if (grad_C.norm() < very_small_value)
    {
        return;
    }

    // Calculate the constraint function value
    const double C = evaluate(V,M);

    // Calculate time-scaled compliance
    const double alpha_tilde = alpha_ / (dt * dt);

    // Inverse mass
    double w0 = 1/(M(indices()[0]));
    double w1 = 1/(M(indices()[1]));
    double w2 = 1/(M(indices()[2]));
    double w3 = 1/(M(indices()[3]));

    const double weighted_sum_of_gradients = (w0 + w1 + w2 + w3)*grad_C.squaredNorm();

    // Calculate \Delta lagrange multiplier
    const double delta_lagrange =
            (-C - alpha_tilde * lagrange) /
            (weighted_sum_of_gradients + alpha_tilde);

    // Update the lagrange multiplier
    lagrange += delta_lagrange;

    // Updates the new values
    V.row(indices()[0]) += w0*delta_lagrange*grad_C.row(0);
    V.row(indices()[1]) += w1*delta_lagrange*grad_C.row(1);
    V.row(indices()[2]) += w2*delta_lagrange*grad_C.row(2);
    V.row(indices()[3]) += w3*delta_lagrange*grad_C.row(3);

}