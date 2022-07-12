//
// Created by aldo on 6/15/22.
//

#include "CXPBDShapeMatchingConstraint.h"

cXPBDShapeMatchingConstraint::scalar_type
cXPBDShapeMatchingConstraint::evaluate(positions_type const& p, masses_type const& M) const
{
    return 0;
}

void cXPBDShapeMatchingConstraint::project(
        positions_type& p,
        positions_type& p0_,
        masses_type const& M,
        scalar_type& lagrange,
        scalar_type const dt,
        gradient_type& F) const
{
    // Calculate the current center of mass
    Eigen::Vector3d x_cm = Eigen::Vector3d::Zero();
    for (int i = 0; i < indices().size(); ++i)
    {
        x_cm += M(this->indices()[i]) * p.row(this->indices()[i]);
    }
    x_cm /= m_total_mass;

    // Calculate A_pq
    Eigen::Matrix3d A_pq = Eigen::Matrix3d::Zero();
    for (int i = 0; i < indices().size(); ++i)
    {
        A_pq += M(this->indices()[i]) * (p.row(this->indices()[i]).transpose() - x_cm) * m_q[i].transpose();
    }

    // Calculate the rotation matrix
    const auto            ATA          = A_pq.transpose() * A_pq;
    const auto            eigen_solver = Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>(ATA);
    const auto            S_inv        = eigen_solver.operatorInverseSqrt();
    const Eigen::Matrix3d R            = A_pq * S_inv;

    assert(R.determinant() > 0);

    // Update the particle positions
    for (int i = 0; i < indices().size(); ++i)
    {
        // Calculate the goal position
        const Eigen::RowVector3d g = R * m_q[i] + x_cm;

        // Move the particle
        p.row(this->indices()[i]) = alpha_ * g + (1.0 - alpha_) * p.row(this->indices()[i]);
    }
}
