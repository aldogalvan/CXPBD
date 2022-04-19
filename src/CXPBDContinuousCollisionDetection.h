//
// Created by aldof on 4/10/2022.
//

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "CXPBDDeformableObject.h"
#include "CXPBDTool.h"

#ifndef CXPBD_CXPBDCONTINUOUSCOLLISIONDETECTION_H
#define CXPBD_CXPBDCONTINUOUSCOLLISIONDETECTION_H

struct ColInfo {
    std::vector<std::set<int>> vertCollisions;
    std::vector<std::set<int>> faceCollisions;

    void resize(size_t s) {
        vertCollisions.resize(s);
        faceCollisions.resize(s);
    }

    void vert_add(size_t off, int index) {
        vertCollisions[off].emplace(index);
    }

    void face_add(size_t off, int index) {
        faceCollisions[off].emplace(index);
    }
};


    ColInfo findCollisions(cXPBDDeformableMesh* def_object,
                           cXPBDTool* tool);

    class CTCD {

    public:

    // Vertex Face continuous collision detection
    // q0: vertex
    // q1-q3: faces
    //
    // For partial credit, the following inputs are tested
    //  1. q0start == q0end, or
    //  2. q1/2/3start == q1/2/3end
    // In other words, either the vertex or the face is stationary.
    static bool vertexFaceCTCD(const Eigen::Vector3d& q0start,
                               const Eigen::Vector3d& q1start,
                               const Eigen::Vector3d& q2start,
                               const Eigen::Vector3d& q3start,
                               const Eigen::Vector3d& q0end,
                               const Eigen::Vector3d& q1end,
                               const Eigen::Vector3d& q2end,
                               const Eigen::Vector3d& q3end,
                               double eta,
                               double& t);

    void Barycentric(Eigen::Vector3d p, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, float& u, float& v, float& w);
    std::tuple<double, double, double, double>
    GetTerms(Eigen::Vector3d& x0, Eigen::Vector3d& y0, Eigen::Vector3d& z0, Eigen::Vector3d& x,
             Eigen::Vector3d& y, Eigen::Vector3d& z);
};




#endif //CXPBD_CXPBDCONTINUOUSCOLLISIONDETECTION_H
