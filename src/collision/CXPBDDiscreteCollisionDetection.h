//
// Created by aldof on 4/28/2022.
//

#include "CXPBDAABB.h"
#include "../world/CXPBDDeformableObject.h"


#ifndef CXPBD_CXPBDDISCRETECOLLISIONDETECTION_H
#define CXPBD_CXPBDDISCRETECOLLISIONDETECTION_H


struct ColInfo {

    int face_index;
    Eigen::Vector3d point;
    set<int> collisionset;

};


//vector<ColInfo*> findCollisions(Eigen::Vector3d p_, double r_ ,cXPBDDeformableMesh& model);

set<int> findCollisions(Eigen::Vector3d p_, double r_ ,cXPBDDeformableMesh& model);

bool sphereTriangleIntersection(Eigen::Vector3d& p_ , double& r_ , Eigen::Vector3d& A, Eigen::Vector3d& B,
                                Eigen::Vector3d& C);

#endif //CXPBD_CXPBDDISCRETECOLLISIONDETECTION_H