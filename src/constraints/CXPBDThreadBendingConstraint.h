//
// Created by agalvan-admin on 8/3/22.
//
#include "CXPBDConstraint.h"

#ifndef CXPBD_CXPBDTHREADBENDINGCONSTRAINT_H
#define CXPBD_CXPBDTHREADBENDINGCONSTRAINT_H

class cXPBDThreadBendingConstraint : public cXPBDConstraint
{
public:
    using self_type      = cXPBDThreadBendingConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDThreadBendingConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta), d_(0.)
    {
        assert(indices.size() == 2u);
        auto const e0 = this->indices()[0];
        auto const e1 = this->indices()[1];

        d_ = (p.row(e0) - p.row(e1)).norm();
    }

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;
    virtual void
    project(positions_type& V, positions_type& V0, masses_type const& M,
            scalar_type& lagrange, scalar_type const dt, gradient_type& F)
    const override;

private:
    scalar_type d_; ///< rest length
};

#endif //CXPBD_CXPBDTHREADBENDINGCONSTRAINT_H
