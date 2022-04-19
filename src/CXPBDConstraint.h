//
// Created by aldof on 4/1/2022.
//

#include "chai3d.h"

#ifndef CXPBD_CXPBDCONSTRAINT_H
#define CXPBD_CXPBDCONSTRAINT_H

using namespace chai3d;

class cXPBDConstraint
{

public:
    using index_type     = std::uint32_t;
    using masses_type    = Eigen::VectorXd;
    using positions_type = Eigen::MatrixXd;
    using gradient_type  = Eigen::Vector3d;
    using scalar_type    = double;

public:
    cXPBDConstraint(std::initializer_list<index_type> indices) : indices_(indices), alpha_(0.0) {}
    cXPBDConstraint(std::initializer_list<index_type> indices, scalar_type const alpha)
            : indices_(indices), alpha_(alpha)
    {
    }

    virtual scalar_type evaluate(positions_type const& V, masses_type const& M) const
    {
        return scalar_type{0.};
    }

    virtual void
    project(positions_type& V, masses_type const& M, scalar_type& lagrange, scalar_type const dt)
    const = 0;
    std::vector<index_type> const& indices() const { return indices_; }

protected:
    scalar_type alpha_;

private:
    std::vector<index_type> indices_;
};
#endif //CXPBD_CXPBDCONSTRAINT_H
