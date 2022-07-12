//
// Created by aldo on 6/15/22.
//

#include "CXPBDConstraint.h"

#ifndef CXPBD_CXPBDSHAPEMATCHINGCONSTRAINT_H
#define CXPBD_CXPBDSHAPEMATCHINGCONSTRAINT_H

class cXPBDShapeMatchingConstraint : public cXPBDConstraint
{
public:
    using self_type      = cXPBDShapeMatchingConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using matrix_type    = std::vector<Eigen::Vector3d>;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDShapeMatchingConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            masses_type const& m,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta), d_(0.)
    {
        auto const e0 = this->indices()[0];
        auto const e1 = this->indices()[1];

        d_ = (p.row(e0) - p.row(e1)).norm();

        // Calculate the initial center of mass and the total mass
        Eigen::Vector3d x_0_cm = Eigen::Vector3d::Zero();
        m_total_mass            = 0.0;
        for (int i = 0; i < indices.size(); ++i)
        {
            m_total_mass += m(this->indices()[i]);
            x_0_cm += m(this->indices()[i]) * p.row(this->indices()[i]);
        }
        x_0_cm /= m_total_mass;

        // Calculate q
        m_q.resize(indices.size());
        for (int i = 0; i < indices.size(); ++i)
        {
            m_q[i] = p.row(this->indices()[i]).transpose() - x_0_cm;
        }
    }

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;
    virtual void
    project(positions_type& V, positions_type& V0, masses_type const& M,
            scalar_type& lagrange, scalar_type const dt, gradient_type& F)
    const override;

private:
    scalar_type d_; ///< rest length
    matrix_type m_q;
    scalar_type m_total_mass;
};

#endif //CXPBD_CXPBDSHAPEMATCHINGCONSTRAINT_H
