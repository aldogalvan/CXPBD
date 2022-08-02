//
// Created by aldo on 6/15/22.
//

#include "CXPBDFixedPointConstraint.h"


 cXPBDFixedPointConstraint::scalar_type cXPBDFixedPointConstraint::evaluate(
         const cXPBDFixedPointConstraint::positions_type &V, const cXPBDFixedPointConstraint::masses_type &M) const {}

void
cXPBDFixedPointConstraint::project(positions_type& V, positions_type& V0, masses_type const& M,
        scalar_type& lagrange, scalar_type const dt, gradient_type& F)  const
{
    V.row(indices()[0]) = p0_;
}