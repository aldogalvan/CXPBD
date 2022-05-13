//
// Created by aldof on 4/10/2022.
//

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "CXPBDDeformableObject.h"
#include "CXPBDTool.h"
#include "CXPBDToolMesh.h"

#ifndef CXPBD_CXPBDCONTINUOUSCOLLISIONDETECTION_H
#define CXPBD_CXPBDCONTINUOUSCOLLISIONDETECTION_H

struct TimeInterval
{
    TimeInterval(double tl, double tu) : l(tl), u(tu)
    {
        if(l > u) std::swap(l, u);
        l = std::max(l, 0.0);
        u = std::min(u, 1.0);
    }

    TimeInterval() : l(0), u(0) {}

    // Returns whether or not the intersection of the intervals is nonempty
    static bool overlap(const TimeInterval &t1, const TimeInterval &t2);
    static bool overlap(const std::vector<TimeInterval> &intervals);

    // Returns the intersection of the intervals **asuming the intersection is nonempty**
    static TimeInterval intersect(const std::vector<TimeInterval> &intervals);

    double l, u;
};


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


    ColInfo findCollisions(Eigen::MatrixXd& p_, Eigen::MatrixXd& plast_, cXPBDDeformableMesh* model,
                           cXPBDToolMesh* tool , double& t , bool& collision);

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


    private:

        // Solves the quadratic equation ax^2 + bx + c = 0, and puts the roots in t0, t1, in ascending order.
        // Returns the number of real roots found.
        static int getQuadRoots(double a, double b, double c, double &t0, double &t1);

        // Conservatively checks if the polynomial of degree degree, with coefficients op, could be positive (if pos is true)
        // or negative (if pos is negative) on the interval [0,1].
        static bool couldHaveRoots(double *op, int degree, bool pos);

        // Looks at the interval [t1, t2], on which a polynomial of degree degree and coefficients op is assumed to have
        // constant sign, and determines if the polynomial is all positive or all negative on that interval.
        // If positive, and pos is true, or if negative, and pos is false, clamps the interval to [0,1] and adds it to
        // intervals.
        static void checkInterval(double t1, double t2, double * op, int degree, std::vector<TimeInterval> &intervals, bool pos);

        // Computes the intervals of x in [0,1] where the polynomial of degree n, with coefficients in op
        // (given in "natural," descending order of power of x) is positive (when pos = true) or negative (if pos = false).
        static void findIntervals(double *op, int n, std::vector<TimeInterval> & intervals, bool pos);

        static void distancePoly3D(const Eigen::Vector3d &x10,
                                   const Eigen::Vector3d &x20,
                                   const Eigen::Vector3d &x30,
                                   const Eigen::Vector3d &v10,
                                   const Eigen::Vector3d &v20,
                                   const Eigen::Vector3d &v30,
                                   double minDSquared,
                                   std::vector<TimeInterval> &result);

        static void barycentricPoly3D(const Eigen::Vector3d &x10,
                                      const Eigen::Vector3d &x20,
                                      const Eigen::Vector3d &x30,
                                      const Eigen::Vector3d &v10,
                                      const Eigen::Vector3d &v20,
                                      const Eigen::Vector3d &v30,
                                      std::vector<TimeInterval> &result);

        static void planePoly3D(const Eigen::Vector3d &x10,
                                const Eigen::Vector3d &x20,
                                const Eigen::Vector3d &x30,
                                const Eigen::Vector3d &v10,
                                const Eigen::Vector3d &v20,
                                const Eigen::Vector3d &v30,
                                std::vector<TimeInterval> &result);

    };



#endif //CXPBD_CXPBDCONTINUOUSCOLLISIONDETECTION_H
