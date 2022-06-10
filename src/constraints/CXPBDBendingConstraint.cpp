//
// Created by aldof on 4/30/2022.
//

#include "CXPBDBendingConstraint.h"


void cXPBDBendingConstraint::project(
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
    auto const v2       = indices.at(2);
    auto const v3       = indices.at(3);
    auto const p0       = p.row(v0);
    auto const p1       = p.row(v1);
    auto const p2       = p.row(v2);
    auto const p3       = p.row(v3);
    auto const w0       = 1. / M(v0);
    auto const w1       = 1. / M(v1);
    auto const w2       = 1. / M(v2);
    auto const w3       = 1. / M(v3);


    gradient_type newcentroid(0, 0, 0);
    for (int j = 0; j < 4u; j++) {
        newcentroid += p.row(indices[j]);
    }

    newcentroid /= 4;
    shape_type B(3, 4);

    for (int j = 0; j < 4u; j++) {
        Eigen::Vector3d newvec = p.row(indices[j]).transpose() - newcentroid;
        B.col(j) = newvec;
    }

    Eigen::Matrix3d M_ = B * A.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(M_, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d R = svd.matrixU() * svd.matrixV().transpose();

    scalar_type const weighted_sum_of_gradients = w0 + w1 + w2 + w3;
    scalar_type const alpha_tilde               = alpha_ / (dt * dt);
    scalar_type const C                         = 1;
    scalar_type const delta_lagrange =
            -(C + alpha_tilde * lagrange) / (weighted_sum_of_gradients + alpha_tilde);

    lagrange += delta_lagrange;
    //p.row(v0) += w0 * R * delta_lagrange;
    //p.row(v1) += w1 * R * delta_lagrange;
    //p.row(v2) += w2 * R * delta_lagrange;
    //p.row(v3) += w3 * R * delta_lagrange;

}