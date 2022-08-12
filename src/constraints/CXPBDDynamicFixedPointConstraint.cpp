//
// Created by agalvan-admin on 8/4/22.
//

#include "CXPBDDynamicFixedPointConstraint.h"

cXPBDDynamicFixedPointConstraint::scalar_type cXPBDDynamicFixedPointConstraint::evaluate(
        const cXPBDDynamicFixedPointConstraint::positions_type &V, const cXPBDDynamicFixedPointConstraint::masses_type &M) const {}

void
cXPBDDynamicFixedPointConstraint::project(positions_type& V, positions_type& V0, masses_type const& M,
                                   scalar_type& lagrange, scalar_type const dt, gradient_type& F)  const
{
    V.row(indices()[0]) = *p0_;
}