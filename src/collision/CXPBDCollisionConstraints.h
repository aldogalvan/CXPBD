
#include <Eigen/Core>
#include "CXPBDDeformableObject.h"

#ifndef CXPBD_CXPBDCOLLISIONCONSTRAINTS_H
#define CXPBD_CXPBDCOLLISIONCONSTRAINTS_H


class cXPBDCollisionConstraints
{
public:

    //! Constructor
    cXPBDCollisionConstraints(cXPBDDeformableMesh* a_object);

    //! Destructor
    ~cXPBDCollisionConstraints();

    //! PUBLIC METHODS
public:

    void proxyCollisionHandling(const Eigen::Vector3d& plast_,  Eigen::Vector3d& p_, const double& t_);

    void implicitCollisionHandling(const Eigen::Vector3d& plast_, const Eigen::Vector3d& p_, const double& t_);

    //! PRIVATE MEMBERS
private:

    cXPBDDeformableMesh* collision_object;

};

#endif //CXPBD_CXPBDCOLLISIONCONSTRAINTS_H