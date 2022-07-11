//
// Created by aldo on 6/10/22.
//

#ifndef CXPBD_CHEROCONSTRAINT_H
#define CXPBD_CHEROCONSTRAINT_H

using namespace chai3d;

class cHEROConstraint
{


public:
    cHEROConstraint(std::initializer_list<uint_32> indices) : indices_(indices), alpha_(0.0), beta_(0.0) {}
    cHEROConstraint(std::initializer_list<uint_32> indices, double const alpha, double const beta)
    : indices_(indices), alpha_(alpha), beta_(beta)
    {
    }

    virtual scalar_type evaluate(cVertexArray const& V, cMassArray const& M) const
    {
        return scalar_type{0.};
    }

    virtual void
    project(cVertexArray& V, cVertexArray& V0, cMassArray const& M,
            scalar_type& lagrange, scalar_type const dt , gradient_type& F)
    const = 0;
    std::vector<index_type> const& indices() const { return indices_; };

    //! PROTECTED MEMBERS

protected:
    scalar_type alpha_;
    scalar_type beta_;
    std::vector<index_type> indices_;
};

#endif //CXPBD_CHEROCONSTRAINT_H
