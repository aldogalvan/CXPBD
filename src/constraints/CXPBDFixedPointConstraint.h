//
// Created by aldo on 6/15/22.
//

#include "CXPBDDeformableObject.h"
#include "CXPBDConstraint.h"

#ifndef CXPBD_CXPBDFIXEDPOINTCONSTRAINT_H
#define CXPBD_CXPBDFIXEDPOINTCONSTRAINT_H


class cXPBDFixedPointConstraint : public cXPBDConstraint
{
public:
    using self_type      = cXPBDFixedPointConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDFixedPointConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta)
    {
        num_vertices = indices.size();
        p0_.resize(num_vertices,3);
        for (int i = 0; i <= num_vertices; i++)
        {
            p0_.row(i) = p.row(this->indices()[i]);
        }

    }

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;
    virtual void
    project(positions_type& V, positions_type& V0, masses_type const& M,
            scalar_type& lagrange, scalar_type const dt, gradient_type& F)
    const override;

private:
    positions_type p0_; ///< rest length
    index_type num_vertices;
};


#endif //CXPBD_CXPBDFIXEDPOINTCONSTRAINT_H
