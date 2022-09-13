#include <Eigen/Dense>
#include "CXPBDConstraint.h"

#ifndef CXPBD_CUDAFIXEDPOINTCONSTRAINT_H
#define CXPBD_CUDAFIXEDPOINTCONSTRAINT_H


class CUDAFixedPointConstraint : public cXPBDConstraint
{
public:
    using self_type      = CUDAFixedPointConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    CUDAFixedPointConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta) , p0_{}
    {
        assert(indices.size() == 1u);
        p0_ = p.row(this->indices()[0]);
    }

    void cuda_project(double* p, double* p0, const double* m,
                      double* lagrange, const double* dt, double* f);

private:

    gradient_type p0_;
};

#endif //CXPBD_CUDAFIXEDPOINTCONSTRAINT_H
