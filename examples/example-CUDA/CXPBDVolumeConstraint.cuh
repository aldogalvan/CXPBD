#include "CXPBDConstraint.h"

class cXPBDVolumeConstraint : public cXPBDConstraint
{
public:
    using self_type      = cXPBDVolumeConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;

public:
    cXPBDVolumeConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0);

    scalar_type evaluate(positions_type const& V, masses_type const& M) const;
    virtual void
    project(positions_type& V, positions_type& Vdot, positions_type& V0, masses_type const& M,
            scalar_type& lagrange, scalar_type const dt)
    const override;

protected:
    scalar_type volume(positions_type const& V) const;

private:
    scalar_type V0_;
};