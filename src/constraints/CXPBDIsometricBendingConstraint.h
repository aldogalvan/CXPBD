//
// Created by aldo on 6/16/22.
//

#include "CXPBDConstraint.h"
#include <Eigen/Dense>

#ifndef CXPBD_CXPBDISOMETRICBENDINGCONSTRAINT_H
#define CXPBD_CXPBDISOMETRICBENDINGCONSTRAINT_H

    class cXPBDIsometricBendingConstraint : public cXPBDConstraint
{


public:
    using self_type      = cXPBDIsometricBendingConstraint;
    using base_type      = cXPBDConstraint;
    using index_type     = std::uint32_t;
    using shape_type     = Eigen::MatrixXd;
    using scalar_type    = double;
    using masses_type    = Eigen::VectorXd;
    using positions_type = typename base_type::positions_type;
    using gradient_type  = typename base_type::gradient_type;
    using hinge_type     = Eigen::MatrixX4d;

public:
        cXPBDIsometricBendingConstraint(
            std::initializer_list<index_type> indices,
            positions_type const& p,
            scalar_type const alpha = 0.0,
            scalar_type const beta = 0.0)
            : base_type(indices, alpha, beta), c_(0.,0.,0.)
    {

        const Eigen::Vector3d& x_0 = p.row(this->indices()[0]);
        const Eigen::Vector3d& x_1 = p.row(this->indices()[1]);
        const Eigen::Vector3d& x_2 = p.row(this->indices()[2]);
        const Eigen::Vector3d& x_3 = p.row(this->indices()[3]);

        const Eigen::Vector3d e0 = x_1 - x_0;
        const Eigen::Vector3d e1 = x_2 - x_1;
        const Eigen::Vector3d e2 = x_0 - x_2;
        const Eigen::Vector3d e3 = x_3 - x_0;
        const Eigen::Vector3d e4 = x_1 - x_3;

        const double cot_01 = calculateCotTheta(e0, -e1);
        const double cot_02 = calculateCotTheta(e0, -e2);
        const double cot_03 = calculateCotTheta(e0, e3);
        const double cot_04 = calculateCotTheta(e0, e4);

        const Eigen::Vector4d K = Eigen::Vector4d(cot_01 + cot_04, cot_02 + cot_03, -cot_01 - cot_02, -cot_03 - cot_04);

        const double A_0 = 0.5 * e0.cross(e1).norm();
        const double A_1 = 0.5 * e0.cross(e3).norm();

        m_Q = (3.0 / (A_0 + A_1)) * K * K.transpose();

    }

    virtual scalar_type evaluate(positions_type const& V, masses_type const& M) const;

    virtual void
    project(positions_type& V, positions_type& V0, masses_type const& M,
            scalar_type& lagrange, scalar_type const dt, gradient_type& F)
    const override;

    inline double calculateCotTheta(const Eigen::Vector3d& x, const Eigen::Vector3d& y);

    inline Eigen::Matrix3d convertVecToCrossOp(const Eigen::Vector3d& vec);



private:
    gradient_type c_;
    shape_type A;
    hinge_type m_Q;
};


#endif //CXPBD_CXPBDISOMETRICBENDINGCONSTRAINT_H
