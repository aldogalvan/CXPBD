//
// Created by aldof on 4/11/2022.
//

/*
#include "chai3d.h"
#include "CXPBDAABB.h"
#include <Eigen/Core>

using namespace chai3d;

#ifndef CXPBD_CXPBDTOOL_H
#define CXPBD_CXPBDTOOL_H

class cXPBDTool : public cToolCursor
{

public:
    using positions_type   = Eigen::MatrixXd;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::MatrixX3d;
    using faces_type       = Eigen::MatrixXi;
    using elements_type    = Eigen::MatrixXi;
    using gradient_type    = Eigen::Vector3d;
    using scalar_type      = double;

    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of cMesh.
    cXPBDTool(cWorld *aParentWorld) : cToolCursor(aParentWorld)
    {

    };

    //! Destructor of cMesh.
    ~cXPBDTool() {};

    //--------------------------------------------------------------------------
    // PUBLIC METHODS - GENERAL
    //--------------------------------------------------------------------------

    gradient_type& positions() { return p_; }
    gradient_type& positions_last() {return plast_;}
    gradient_type& velocity() { return v_; }
    scalar_type const & mass() { return m_; }
    scalar_type const & radius() {return r_;}

    // This method connects the mesh to chai3d
    void connectToChai3d(void);

    // This method updates the chai3d mesh
    void updateChai3d(void);

    // This method translates the object
    void setLocalPos(Eigen::Vector3d a_pos);

    // This method scales the object
    void scaleObject(double a_scale);

    // This method sets the velocities of the object
    void setVelocities(velocities_type a_velocity)
    {
        v_ = a_velocity;
    };


    // This method sets the mass of the deformable object
    void setMass(double a_mass)
    {
        m_ = a_mass;
    }


    void setHapticPos(Eigen::Vector3d a_pos)
    {
        last_proxy_pos = proxy_pos;
        proxy_pos = a_pos;
    }

    void updatePos()
    {
        plast_ = p_;
        setLocalPos(proxy_pos);
    }


    void projectPos(double t)
    {
        p_ = plast_ + t*(p_ - plast_);
    }

public:
    bool collision;
    gradient_type normal;
    double depth;

private:
    gradient_type p_;            ///< Rest positions
    gradient_type plast_;         ///< Last positions
    gradient_type v_;            ///< Per-vertex velocity
    gradient_type vlast_;        ///< Last velocity
    gradient_type proxy_pos;      ///< Haptic position
    gradient_type last_proxy_pos; ///< Last haptic position
    scalar_type m_;                ///< Per-vertex mass_when_unfixed
    scalar_type r_;                ///< Radius of the sphere
};

#endif //CXPBD_CXPBDTOOL_H
*/