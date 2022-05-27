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
        positions_type& v,
        positions_type& p0_,
        masses_type const& M,
        scalar_type& lagrange,
        scalar_type const dt) const
{
    auto const& indices = this->indices();
    auto const v0       = indices.at(0);
    auto const v1       = indices.at(1);
    auto const p0       = p.row(v0);
    auto const p1       = p.row(v1);
    auto const pdot0    = v.row(v0);
    auto const pdot1    = v.row(v1);
    auto const w0       = 1. / M(v0);
    auto const w1       = 1. / M(v1);
    auto const n        = (p0 - p1).normalized();
    auto const C        = evaluate(p, M);


    // <n,n> = 1 and <-n,-n> = 1
    scalar_type const weighted_sum_of_gradients = w0 + w1;
    scalar_type const alpha_tilde               = alpha_ / (dt * dt);
    scalar_type const beta_tilde                = beta_ * dt * dt;
    scalar_type const gamma                     = alpha_tilde * beta_tilde / dt;
    scalar_type const delta_lagrange_0 =
            (-C - alpha_tilde * lagrange - gamma * n.dot(pdot0) * dt) /
            ((1 + gamma)*weighted_sum_of_gradients + alpha_tilde);
    scalar_type const delta_lagrange_1 =
            (-C - alpha_tilde * lagrange - gamma * -n.dot(pdot0) * dt) /
            ((1 + gamma)*weighted_sum_of_gradients + alpha_tilde);

    lagrange += delta_lagrange_0;
    lagrange += delta_lagrange_1;

    p.row(v0) += w0 * n * delta_lagrange_0;
    p.row(v1) += w1 * -n * delta_lagrange_1;
}
