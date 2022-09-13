#include <Eigen/Dense>
#include "CXPBDConstraint.h"

#ifndef CXPBD_CUDAVOLUMECONSTRAINT_H
#define CXPBD_CUDAVOLUMECONSTRAINT_H

#endif //CXPBD_CUDAVOLUMECONSTRAINT_H


class CUDAVolumeConstraint : public cXPBDConstraint
{
public:
    using self_type      = CUDAVolumeConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:

    CUDAVolumeConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0);

    void cuda_project(double* p, double* p0, const double* m,
                      double* lagrange, const double* dt, double* f);

protected:
    scalar_type volume(positions_type const& V) const;

public:
    scalar_type V0_;
};
