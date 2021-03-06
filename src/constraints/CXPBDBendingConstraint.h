//
// Created by aldof on 4/30/2022.
//

#include "CXPBDConstraint.h"

#ifndef CXPBD_CXPBDBENDINGCONSTRAINT_H
#define CXPBD_CXPBDBENDINGCONSTRAINT_H

class cXPBDBendingConstraint : public cXPBDConstraint
{


public:
    using self_type      = cXPBDBendingConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using shape_type     = Eigen::MatrixXd;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDBendingConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const di_ang = 0.0,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta), m_di_ang(di_ang) , c_(0.,0.,0.)
    {
        
        assert(indices.size() == 4u);
        auto const h0 = this->indices()[0];
        auto const h1 = this->indices()[1];
        auto const h2 = this->indices()[2];
        auto const h3 = this->indices()[3];

        c_ = (p.row(h0) + p.row(h1) + p.row(h2) + p.row(h3)) / 4;

        A.resize(3,4);

        for (int j = 0; j < 4u; j++) {
            Eigen::Vector3d oldvec = p.row(this->indices()[j]).transpose() - c_;
            A.col(j) = oldvec;
        }

    }

    virtual scalar_type evaluate(positions_type const& V, masses_type const& M) const;

    virtual void
    project(positions_type& V, positions_type& V0, masses_type const& M,
            scalar_type& lagrange, scalar_type const dt, gradient_type& F)
    const override;


private:
    gradient_type c_;
    shape_type A;
    scalar_type m_di_ang;
};


#endif //CXPBD_CXPBDBENDINGCONSTRAINT_H
