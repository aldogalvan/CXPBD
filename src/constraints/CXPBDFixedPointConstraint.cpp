//
// Created by aldo on 6/15/22.
//

#include "CXPBDFixedPointConstraint.h"

void
cXPBDFixedPointConstraint::project(positions_type& V, positions_type& V0, masses_type const& M,
        scalar_type& lagrange, scalar_type const dt, gradient_type& F)  const
{

    for (int i = 0u; i < num_vertices ; i++)
    {
        V.row(i) = p0_.row(indices()[i]);
    }
}