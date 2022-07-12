//
// Created by aldof on 4/28/2022.
//

#include "CXPBDDiscreteCollisionDetection.h"

using namespace std;
using namespace Eigen;

/*
bool sphereTriangleIntersection(Eigen::Vector3d& p_ , double& r_ , Eigen::Vector3d& A, Eigen::Vector3d& B,
                                Eigen::Vector3d& C , Eigen::Vector3d n_, ColInfo& collisions)
{
    A -= p_;
    B -= p_;
    C -= p_;

    //First check if sphere is close to triangle plane

    double d = A.dot(n_);
    double e = n_.dot(n_);
    double rr = r_ * r_;

    // If so then no collision
    if (d*d > rr *e)
        return false;

    double aa = A.squaredNorm();
    double ab = A.dot(B);
    double ac = A.dot(C);
    double bb = B.squaredNorm();
    double bc = B.dot(C);
    double cc = C.squaredNorm();


}
 */
static void Barycentric(float out[3], const Vector3d& A, const Vector3d& B, const Vector3d& Q)
{
    Vector3d AB = B - A;
    Vector3d QA = A - Q;
    Vector3d QB = B - Q;

    float divisor = AB.dot(AB);

    out[0] = QB.dot(AB);
    out[1] = -QA.dot(AB);
    out[2] = divisor;
}

// Convert a point Q from Cartesian coordinates to Barycentric coordinates (u, v, w)
// with respect to a triangle ABC.
// The last output value is the divisor.
static void Barycentric(float out[4], const Vector3d& A, const Vector3d& B, const Vector3d& C,
                        const Vector3d& Q)
{
    Vector3d AB = B - A;
    Vector3d AC = C - A;

    Vector3d QA = A - Q;
    Vector3d QB = B - Q;
    Vector3d QC = C - Q;

    Vector3d QB_x_QC = QB.cross(QC);
    Vector3d QC_x_QA = QC.cross(QA);
    Vector3d QA_x_QB = QA.cross(QB);

    Vector3d AB_x_AC = AB.cross(AC);

    float divisor = AB_x_AC.dot(AB_x_AC);

    out[0] = QB_x_QC.dot(AB_x_AC);
    out[1] = QC_x_QA.dot(AB_x_AC);
    out[2] = QA_x_QB.dot(AB_x_AC);
    out[3] = divisor;
}


// Return the closest point on a triangle ABC to a point Q.
bool SolvePoint(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& Q,
           Vector3d& p)
{
    // Test vertex regions
    float wAB[3], wBC[3], wCA[3];
    Barycentric(wAB, A, B, Q);
    Barycentric(wBC, B, C, Q);
    Barycentric(wCA, C, A, Q);

    // R A
    if (wAB[1] <= 0.0f && wCA[0] <= 0.0f)
    {
        p = A;
        return 1;
    }

    // R B
    if (wAB[0] <= 0.0f && wBC[1] <= 0.0f)
    {
        p = B;
        return 1;
    }

    // R C
    if (wBC[0] <= 0.0f && wCA[1] <= 0.0f)
    {
        p = C;
        return 1;
    }

    // Test edge regions
    float wABC[4];
    Barycentric(wABC, A, B, C, Q);

    // This is used to help testing if the face degenerates
    // into an edge.
    float area = wABC[3];

    // R AB
    if (wAB[0] > 0.0f && wAB[1] > 0.0f && area * wABC[2] <= 0.0f)
    {
        float s = wAB[2];
        assert(s > 0.0f);
        p = (wAB[0] * A + wAB[1] * B) / s;
        return 1;
    }

    // R BC
    if (wBC[0] > 0.0f && wBC[1] > 0.0f && area * wABC[0] <= 0.0f)
    {
        float s = wBC[2];
        assert(s > 0.0f);
        p = (wBC[0] * B + wBC[1] * C) / s;
        return 1;
    }

    // R CA
    if (wCA[0] > 0.0f && wCA[1] > 0.0f && area * wABC[1] <= 0.0f)
    {
        float s = wCA[2];
        assert(s > 0.0f);
        p = (wCA[0] * C + wCA[1] * A) / s;
        return 1;
    }

    // R ABC/ACB
    float s = wABC[3];
    if (s <= 0.0f)
    {
        // Give up. Triangle is degenerate.
        // You should handle all cases correctly in production code.
        return 0;
    }

    assert(s > 0.0f);
    p = (wABC[0] * A + wABC[1] * B + wABC[2] * C) / s;
    return 1;
}

set<int> findCollisions(Eigen::Vector3d p_, double r_, cXPBDDeformableMesh& model)
{

    set<int> ret;

    auto& pdes = model.positions();
    model.buildAABBBoundaryBox();
    auto toolBB_ = buildAABB(p_,r_);
    const auto& f_ = model.faces();
    const auto& N_ = model.normals();
    std::vector<Collision> potentialCollisions;
    intersect(toolBB_, model.bb(), potentialCollisions);

    for (const auto& it : potentialCollisions)
    {

        Eigen::Vector3i face = f_.row(it.collidingTriangle2);
        Eigen::Vector3d A = pdes.row(face(0));
        Eigen::Vector3d B = pdes.row(face(1));
        Eigen::Vector3d C = pdes.row(face(2));
        Eigen::Vector3d normal = N_.row(it.collidingTriangle2);
        Eigen::Vector3d p;

        if ( SolvePoint(A,B,C,p_,p))
        {
            if ((p_-p).squaredNorm() < r_*r_)
            {

                //ColInfo* newcollision;
                //newcollision  = new ColInfo;
                //std::cout << "GOOD" << std::endl;
                //newcollision->face_index = it.collidingTriangle2;
                //newcollision->point = p;
                ret.insert(face(0));
                ret.insert(face(1));
                ret.insert(face(2));
                //std::cout << "GOOD" << std::endl;
            }
        }
    }
    return ret;
}
