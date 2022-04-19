//
// Created by aldof on 4/2/2022.
//

#ifndef CXPBD_CXPBDNEOHOOKEANCONSTRAINT_H
#define CXPBD_CXPBDNEOHOOKEANCONSTRAINT_H

#include "CXPBDConstraint.h"


class cXPBDNeoHookeanConstraint : public cXPBDConstraint {
public:
    using self_type = cXPBDNeoHookeanConstraint;
    using base_type = cXPBDConstraint;
    using index_type = std::uint32_t;
    using scalar_type = typename cXPBDConstraint::scalar_type;
    using masses_type = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type = typename base_type::gradient_type;

public:
    cXPBDNeoHookeanConstraint(
            std::initializer_list<index_type> indices,
            positions_type const &p,
            scalar_type young_modulus,
            scalar_type poisson_ratio,
            scalar_type const alpha = 0.0);

    virtual void
    project(positions_type &p, masses_type const &m, scalar_type &lagrange, scalar_type const dt)
    const override;

protected:
    scalar_type signed_volume(positions_type const &V) const;

private:
    scalar_type V0_;
    Eigen::Matrix3d DmInv_;
    scalar_type mu_;
    scalar_type lambda_;
};

#endif //CXPBD_CXPBDNEOHOOKEANCONSTRAINT_H
