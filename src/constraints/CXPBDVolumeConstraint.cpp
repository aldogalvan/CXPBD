//
// Created by aldof on 4/2/2022.
//

#include "CXPBDVolumeConstraint.h"

cXPBDVolumeConstraint::cXPBDVolumeConstraint(
        std::initializer_list<index_type> indices,
        positions_type const& p,
        scalar_type const alpha,
        scalar_type const beta)
        : base_type(indices, alpha, beta), V0_{0.}
{
    assert(indices.size() == 4u);
    V0_ = volume(p);
}

cXPBDVolumeConstraint::scalar_type
cXPBDVolumeConstraint::volume(positions_type const& V) const
{
    Eigen::RowVector3d const p0 = V.row(indices()[0]);
    Eigen::RowVector3d const p1 = V.row(indices()[1]);
    Eigen::RowVector3d const p2 = V.row(indices()[2]);
    Eigen::RowVector3d const p3 = V.row(indices()[3]);

    auto const vol = (1. / 6.) * (p1 - p0).cross(p2 - p0).dot(p3 - p0);
    return std::abs(vol);
}

cXPBDVolumeConstraint::scalar_type
cXPBDVolumeConstraint::evaluate(positions_type const& p, masses_type const& m) const
{
    return volume(p) - V0_;
}

void cXPBDVolumeConstraint::project(
        positions_type& p,
        positions_type& p0_,
        masses_type const& m,
        scalar_type& lagrange,
        scalar_type const dt,
        gradient_type& F) const
{
    auto const v0 = indices()[0];
    auto const v1 = indices()[1];
    auto const v2 = indices()[2];
    auto const v3 = indices()[3];

    auto const w0 = 1. / m(v0);
    auto const w1 = 1. / m(v1);
    auto const w2 = 1. / m(v2);
    auto const w3 = 1. / m(v3);

    Eigen::RowVector3d const p0 = p.row(v0);
    Eigen::RowVector3d const p1 = p.row(v1);
    Eigen::RowVector3d const p2 = p.row(v2);
    Eigen::RowVector3d const p3 = p.row(v3);

    Eigen::Vector3d const pdot0  = p0 - p0_.row(v0);
    Eigen::Vector3d const pdot1  = p1 - p0_.row(v1);
    Eigen::Vector3d const pdot2  = p2 - p0_.row(v2);
    Eigen::Vector3d const pdot3  = p3 - p0_.row(v3);

    auto const C = evaluate(p, m);

    Eigen::RowVector3d const grad0 = (1. / 6.) * (p1 - p2).cross(p3 - p2);
    Eigen::RowVector3d const grad1 = (1. / 6.) * (p2 - p0).cross(p3 - p0);
    Eigen::RowVector3d const grad2 = (1. / 6.) * (p0 - p1).cross(p3 - p1);
    Eigen::RowVector3d const grad3 = (1. / 6.) * (p1 - p0).cross(p2 - p0);

    auto const weighted_sum_of_gradients = w0 * grad0.squaredNorm() + w1 * grad1.squaredNorm() +
                                           w2 * grad2.squaredNorm() + w3 * grad3.squaredNorm();

    scalar_type const alpha_tilde = alpha_ / (dt * dt);
    scalar_type const beta_tilde = beta_ * dt * dt;
    scalar_type const gamma = (alpha_tilde * beta_tilde ) / dt;

    if (abs(alpha_tilde+weighted_sum_of_gradients) < 1e-6)
        return;

    scalar_type const delta_lagrange =
            (-C - alpha_tilde * lagrange - gamma * (grad0.dot(pdot0) + grad1.dot(pdot1)
            + grad2.dot(pdot2) + grad3.dot(pdot3))) / ((1 + gamma)*weighted_sum_of_gradients + alpha_tilde);

    lagrange += delta_lagrange;
    F += delta_lagrange*(grad0 + grad1 + grad2 + grad3);


    p.row(v0) += w0 * grad0 * delta_lagrange;
    p.row(v1) += w1 * grad1 * delta_lagrange;
    p.row(v2) += w2 * grad2 * delta_lagrange;
    p.row(v3) += w3 * grad3 * delta_lagrange;

}