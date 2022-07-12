//
// Created by aldof on 4/2/2022.
//

#include "CXPBDEdgeLengthConstraint.h"

cXPBDEdgeConstraint::scalar_type
cXPBDEdgeConstraint::evaluate(positions_type const& p, masses_type const& M) const
{
    auto const v0 = indices().at(0);
    auto const v1 = indices().at(1);
    auto const p0 = p.row(v0);
    auto const p1 = p.row(v1);

    return (p0 - p1).norm() - d_;
}

void cXPBDEdgeConstraint::project(
        positions_type& p,
        positions_type& p0_,
        masses_type const& M,
        scalar_type& lagrange,
        scalar_type const dt,
        gradient_type& F) const
{
    auto const& indices = this->indices();
    auto const v0       = indices.at(0);
    auto const v1       = indices.at(1);
    auto const p0         = p.row(v0);
    auto const p1         = p.row(v1);
    auto const pdot0     = p0 - p0_.row(v0);
    auto const pdot1     = p1 - p0_.row(v1);
    auto const w0          = 1. / M(v0);
    auto const w1          = 1. / M(v1);
    Eigen::Vector3d n = (p0 - p1);
    double d = n.norm();
    if (d > 1e-6)
        n = n / d;
    else
        return;

    auto const C        = evaluate(p, M);


    // <n,n> = 1 and <-n,-n> = 1
    auto const grad0             = n;
    auto const grad1         = -n;
    scalar_type const weighted_sum_of_gradients = w0 + w1;
    scalar_type const alpha_tilde               = alpha_ / (dt * dt);
    scalar_type const beta_tilde                = beta_ * dt * dt;
    scalar_type const gamma                     = alpha_tilde * beta_tilde / dt;


    if (alpha_tilde + weighted_sum_of_gradients < 1e-6)
        return;


    scalar_type const delta_lagrange =
            (-C - alpha_tilde * lagrange - gamma * (grad0.dot(pdot0) + grad1.dot(pdot1))) /
            ((1 + gamma)*weighted_sum_of_gradients + alpha_tilde);

    lagrange += delta_lagrange;

    F += delta_lagrange*(grad0 + grad1);

    p.row(v0) += w0 * grad0 * delta_lagrange;
    p.row(v1) += w1 * grad1 * delta_lagrange;
}
