#include <Eigen/Dense>
#include "CXPBDConstraint.h"

#ifndef CXPBD_CUDAEDGECONSTRAINT_H
#define CXPBD_CUDAEDGECONSTRAINT_H

class CUDAEdgeConstraint : public cXPBDConstraint
{
public:
    using self_type      = CUDAEdgeConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    CUDAEdgeConstraint(
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

    void cuda_project(double* p, double* p0, const double* m,
            double* lagrange, const double* dt, double* f);

private:
    scalar_type d_; ///< rest length
};

#endif //CXPBD_CUDAEDGECONSTRAINT_H
