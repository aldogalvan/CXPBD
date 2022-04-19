#ifndef PSIM_CORE_SHARED_VECTOR_MATH_H
#define PSIM_CORE_SHARED_VECTOR_MATH_H

#include <Eigen/Core>

class VectorMath
{
public:
    static const Eigen::Matrix3d crossProductMatrix(const Eigen::Vector3d &v);
    static const Eigen::Matrix3d rotationMatrix(const Eigen::Vector3d &axisAngle);
    static const Eigen::Vector3d axisAngle(const Eigen::Matrix3d &rotationMatrix);
    static const Eigen::Vector3d perpToAxis(const Eigen::Vector3d &v);
    static const Eigen::Matrix3d DrotVector(const Eigen::Vector3d &axisangle, const Eigen::Vector3d &rotatingVector);
    static const Eigen::Matrix3d TMatrix(const Eigen::Vector3d &v);
    static double randomUnitIntervalReal();
    static const Eigen::Vector3d randomPointOnSphere();
};

#ifndef M_PI
// This fixes Windows build, but be aware EIGEN_PI is not part of its API
// and subject to change.
//
// A long term solution is to use coming C++ 20's std::numbers::pi
constexpr double M_PI = EIGEN_PI;
#endif

#endif // MATH_H
