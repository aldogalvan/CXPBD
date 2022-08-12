//
// Created by agalvan-admin on 8/4/22.
//

#include "CXPBDDeformableObject.h"
#include "CXPBDConstraint.h"

#ifndef CXPBD_CXPBDDYNAMICFIXEDPOINTCONSTRAINT_H
#define CXPBD_CXPBDDYNAMICFIXEDPOINTCONSTRAINT_H

class cXPBDDynamicFixedPointConstraint : public cXPBDConstraint
{
public:
    using self_type      = cXPBDDynamicFixedPointConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDDynamicFixedPointConstraint(
            std::initializer_list<index_type> indices,
            shared_ptr<gradient_type> p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta) , p0_{}
    {
        assert(indices.size() == 1u);
        p0_ = p;
    }

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;

    virtual void project(positions_type& V, positions_type& V0, masses_type const& M,
                         scalar_type& lagrange, scalar_type const dt, gradient_type& F)
    const override;

private:

    shared_ptr<gradient_type> p0_;
};

#endif //CXPBD_CXPBDDYNAMICFIXEDPOINTCONSTRAINT_H
