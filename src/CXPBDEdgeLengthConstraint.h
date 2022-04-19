//
// Created by aldof on 4/2/2022.
//

#include "CXPBDConstraint.h"

#ifndef CXPBD_CXPBDEDGELENGTHCONSTRAINT_H
#define CXPBD_CXPBDEDGELENGTHCONSTRAINT_H

class cXPBDEdgeConstraint : public cXPBDConstraint
{
public:
    using self_type      = cXPBDEdgeConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDEdgeConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0)
            : base_type(indices, alpha), d_(0.)
    {
        assert(indices.size() == 2u);
        auto const e0 = this->indices()[0];
        auto const e1 = this->indices()[1];

        d_ = (p.row(e0) - p.row(e1)).norm();
    }

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;
    virtual void
    project(positions_type& V, masses_type const& M, scalar_type& lagrange, scalar_type const dt)
    const override;

private:
    scalar_type d_; ///< rest length
};

#endif //CXPBD_CXPBDEDGELENGTHCONSTRAINT_H
