//
// Created by aldof on 4/10/2022.
//

#include <memory>
#include "CXPBDContinuousCollisionDetection.h"
#include "../shared/rpoly.h"
#include "CXPBDAABB.h"

using namespace std;

bool TimeInterval::overlap(const TimeInterval &t1, const TimeInterval &t2)
{
    return !(t1.l > t2.u || t2.l > t1.u);
}

bool TimeInterval::overlap(const std::vector<TimeInterval> &intervals)
{
    for(std::vector<TimeInterval>::const_iterator it1 = intervals.begin(); it1 != intervals.end(); ++it1)
    {
        std::vector<TimeInterval>::const_iterator it2 = it1;
        for(++it2; it2 != intervals.end(); ++it2)
            if(!overlap(*it1, *it2))
                return false;
    }
    return true;
}

TimeInterval TimeInterval::intersect(const std::vector<TimeInterval> &intervals)
{
    TimeInterval isect(0.0, 1.0);
    for(std::vector<TimeInterval>::const_iterator it = intervals.begin(); it != intervals.end(); ++it)
    {
        isect.l = max(it->l, isect.l);
        isect.u = min(it->u, isect.u);
    }
    return isect;
}

int CTCD::getQuadRoots(double a, double b, double c, double &t0, double &t1) {
    int roots = 0;
    int sign = 1;

    if (b < 0)
        sign = -1;

    double D = b * b - 4 * a * c;
    if (D >= 0)
    {
        roots = 2;
        double q = -0.5 * (b + sign * sqrt(D));
        t0 = q / a;
        t1 = c / q;
        if (t0 > t1)
            std::swap(t1, t0);
    }
    return roots;
}


void CTCD::checkInterval(double t1, double t2, double * op, int degree, vector<TimeInterval> &intervals, bool pos)
{
    // clamp values
    t1 = max(0.0, t1);
    t2 = max(0.0, t2);
    t1 = min(1.0, t1);
    t2 = min(1.0, t2);

    double tmid = (t2 + t1) / 2;
    double f = op[0];
    for(int i=1; i<=degree; i++)
    {
        f *= tmid;
        f += op[i];
    }

    if (pos && f >= 0)
        intervals.push_back(TimeInterval(t1, t2));
    else if (!pos && f <= 0)
        intervals.push_back(TimeInterval(t1, t2));
}

bool CTCD::couldHaveRoots(double *op, int degree, bool pos) {
    double result = 0;
    if ((pos && op[0] > 0) || (!pos && op[0] < 0))
        result = op[0];
    for (int i = 1; i < degree; i++)
    {
        result *= 1.0;
        if ((pos && op[i] > 0) || (!pos && op[i] < 0))
            result += op[i];
    }
    result *= 1.0;
    result += op[degree];
    return !((pos && result < 0) || (!pos && result > 0));
}



void CTCD::findIntervals(double *op, int n, vector<TimeInterval> & intervals, bool pos)
{
    int roots=0;
    int reducedDegree=n;

    if(n>6)
    {
        assert(!"Polynomials of degree > 6 not supported");
        return;
    }

    double time[6];

    // We don't care one bit about these imaginary roots
    double zeroi[6];

    // normalize
    double maxval = 0;
    for(int i=0; i<=n; i++)
        maxval = std::max(maxval, fabs(op[i]));
    if(maxval != 0)
        for(int i=0; i<=n; i++)
            op[i] /= maxval;

    for (int i = 0; i < n; i++)
    {
        if (op[i] == 0)
            reducedDegree--;
        else
            break;
    }

    if (reducedDegree < n)
    {
        for (int i = 0; i <= reducedDegree; i++)
            op[i] = op[i + n - reducedDegree];
    }

    if (reducedDegree > 2) {
        if (!couldHaveRoots(op, reducedDegree, pos))
            return;

        RootFinder rf;

        roots = rf.rpoly(op, reducedDegree, time, zeroi);
    }
    else if (reducedDegree == 2)
    {
        roots = getQuadRoots(op[0], op[1], op[2], time[0], time[1]);
    }
    else if (reducedDegree == 1)
    {
        time[0] = -op[1] / op[0];
        roots = 1;
    }
    else
    {
        // both points stationary -- check if colliding at t=0
        if ((!pos && op[0] <= 0) || (pos && op[0] >= 0))
            intervals.push_back(TimeInterval(0, 1.0));
        return;
    }

    // check intervals
    if (roots > 0)
    {
        std::sort(time, time + roots);
        if (time[0] >= 0)
            checkInterval(0, time[0], op, reducedDegree, intervals, pos);
        for (int i = 0; i < roots - 1; i++) {
            if (!((time[i] < 0 && time[i + 1] < 0) || (time[i] > 1.0 && time[i + 1] > 1.0)))
                checkInterval(time[i], time[i + 1], op, reducedDegree, intervals, pos);
        }
        if (time[roots - 1] <= 1.0)
            checkInterval(time[roots - 1], 1.0, op, reducedDegree, intervals, pos);
    }
    else
    {
        checkInterval(0.0, 1.0, op, reducedDegree, intervals, pos);
    }
}

void CTCD::barycentricPoly3D(const Eigen::Vector3d &x10,
                             const Eigen::Vector3d &x20,
                             const Eigen::Vector3d &x30,
                             const Eigen::Vector3d &v10,
                             const Eigen::Vector3d &v20,
                             const Eigen::Vector3d &v30,
                             vector<TimeInterval> &result)
{
    // alpha > 0
    double A = x10.dot(x10);
    double B = 2 * x10.dot(v10);
    double C = (v10).dot(v10);
    // e0.e1
    double D = (x20).dot(x10);
    double E = (x20).dot(v10) + (v20).dot(x10);
    double F = (v20).dot(v10);
    //(q0-q1).e0
    double G = (x30).dot(x20);
    double H = (x30).dot(v20) + (v30).dot(x20);
    double I = (v30).dot(v20);
    //(q0-q1).e1
    double J = (x30).dot(x10);
    double K = (x30).dot(v10) + (v30).dot(x10);
    double L = (v30).dot(v10);

    double op[5];

    op[0] = F * L - C * I;
    op[1] = F * K + E * L - C * H - B * I;
    op[2] = F * J + D * L + E * K - C * G - A * I - B * H;
    op[3] = D * K + E * J - A * H - B * G;
    op[4] = D * J - A * G;

    findIntervals(op, 4, result, true);
}

void CTCD::planePoly3D(const Eigen::Vector3d &x10,
                       const Eigen::Vector3d &x20,
                       const Eigen::Vector3d &x30,
                       const Eigen::Vector3d &v10,
                       const Eigen::Vector3d &v20,
                       const Eigen::Vector3d &v30,
                       vector<TimeInterval> &result)
{
    double op[4];
    op[0] = v10.dot(v20.cross(v30));
    op[1] = x10.dot(v20.cross(v30)) + v10.dot(x20.cross(v30)) + v10.dot(v20.cross(x30));
    op[2] = x10.dot(x20.cross(v30)) + x10.dot(v20.cross(x30)) + v10.dot(x20.cross(x30));
    op[3] = x10.dot(x20.cross(x30));
    findIntervals(op, 3, result, true);
}

void CTCD::distancePoly3D(const Eigen::Vector3d &x10,
                          const Eigen::Vector3d &x20,
                          const Eigen::Vector3d &x30,
                          const Eigen::Vector3d &v10,
                          const Eigen::Vector3d &v20,
                          const Eigen::Vector3d &v30,
                          double minDSquared,
                          vector<TimeInterval> &result)
{
    double A = v10.dot(v20.cross(v30));
    double B = x10.dot(v20.cross(v30)) + v10.dot(x20.cross(v30)) + v10.dot(v20.cross(x30));
    double C = x10.dot(x20.cross(v30)) + x10.dot(v20.cross(x30)) + v10.dot(x20.cross(x30));
    double D = x10.dot(x20.cross(x30));
    Eigen::Vector3d E = x20.cross(x30);
    Eigen::Vector3d F = x20.cross(v30) + v20.cross(x30);
    Eigen::Vector3d G = v20.cross(v30);

    double op[7];
    op[0] = A * A;
    op[1] = 2 * A * B;
    op[2] = B * B + 2 * A * C - G.dot(G) * minDSquared;
    op[3] = 2 * A * D + 2 * B * C - 2 * G.dot(F) * minDSquared;
    op[4] = 2 * B * D + C * C - (2 * G.dot(E) + F.dot(F)) * minDSquared;
    op[5] = 2 * C * D - 2 * F.dot(E) * minDSquared;
    op[6] = D * D - E.dot(E) * minDSquared;
    findIntervals(op, 6, result, false);
}

ColInfo findCollisions(Eigen::MatrixXd& pgoal_, Eigen::Vector3d& p_, Eigen::Vector3d& plast_, cXPBDDeformableMesh& model, double &t_)
{
    ColInfo ret;

    auto& faceCollisions = ret.faceCollisions;

    // ONLY HAVE ONE PONT REPRESENTING THE DEVICE
    // vertCollisions.resize(1);
    // faceCollisions.resize(model.numFaces());
    double radius = 0.01;

    model.buildAABBBoundaryBox(pgoal_);
    auto toolBB_ = buildAABB(plast_,p_,radius);

    const auto& f_ = model.faces();
    const auto& v_ = model.positions();
    const auto& vlast_ = model.positions();

    std::vector<Collision> potentialCollisions;
    intersect(toolBB_, model.bb(), potentialCollisions);
    int it_count = 0;

    for (const auto& it : potentialCollisions) {

        Eigen::Vector3i face = f_.row(it.collidingTriangle2);
        Eigen::Vector3d alast_ = vlast_.row(face(0));
        Eigen::Vector3d blast_ = vlast_.row(face(1));
        Eigen::Vector3d clast_  = vlast_.row(face(2));
        Eigen::Vector3d a_ = v_.row(face(0));
        Eigen::Vector3d b_ = v_.row(face(1));
        Eigen::Vector3d c_ = v_.row(face(2));
        Eigen::Vector3d vlastnormal_ = (blast_ - alast_).cross(clast_-alast_).normalized();
        Eigen::Vector3d vnormal_ = (b_ - a_).cross(c_ - a_).normalized();


        if (CTCD::vertexFaceCTCD(
                (plast_ ), // - vlastnormal_*radius
                alast_.transpose(),
                blast_.transpose(),
                clast_.transpose(),
                (p_ ), //- vnormal_*radius
                a_.transpose(),
                b_.transpose(),
                c_.transpose(),
                1e-6,
                t_))
        {

            ret.collision = true;
            faceCollisions.insert(it.collidingTriangle2);
            //ret.col.insert(std::make_tuple(it.collidingTriangle2,vlastnormal_,vnormal_));
            it_count++;
            return ret;
        }

    }

    /*
    for (const auto& it : potentialCollisions) {
        Eigen::Vector3i face = F_tool.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = F_(it.collidingTriangle2, vert);
            double t;

            if (CTCD::vertexFaceCTCD(
                    plast_.row(vidx).transpose(),
                    V_tool_last.row(face(0)).transpose(),
                    V_tool_last.row(face(1)).transpose(),
                    V_tool_last.row(face(2)).transpose(),
                    p_.row(vidx).transpose(),
                    V_tool.row(face(0)).transpose(),
                    V_tool.row(face(1)).transpose(),
                    V_tool.row(face(2)).transpose(),
                    1e-6,
                    t))
            {
                faceCollisions[it.collidingTriangle1].emplace(vidx);
            }
        }
    }
     */

    return ret;

}


// TODO: Implmentation of Vertex-Face CTCD Narrow Phase
bool CTCD::vertexFaceCTCD(const Eigen::Vector3d& q0start,
                          const Eigen::Vector3d& q1start,
                          const Eigen::Vector3d& q2start,
                          const Eigen::Vector3d& q3start,
                          const Eigen::Vector3d& q0end,
                          const Eigen::Vector3d& q1end,
                          const Eigen::Vector3d& q2end,
                          const Eigen::Vector3d& q3end,
                          double eta,
                          double& t)
{

    double mind = eta*eta;

    Eigen::Vector3d v0 = q0end - q0start;
    Eigen::Vector3d v1 = q1end - q1start;
    Eigen::Vector3d v2 = q2end - q2start;
    Eigen::Vector3d v3 = q3end - q3start;

    // time intervals during which v is colinear with the edge,
    // on the side of e1 towards e2, and on the side of e2 towards e1
    vector<TimeInterval> coplane, e1, e2, e3;

    Eigen::Vector3d x10 = q0start - q1start;
    Eigen::Vector3d v10 = v0 - v1;
    Eigen::Vector3d x20 = (q3start - q1start).cross(q2start - q1start);
    Eigen::Vector3d v20 = (v3 - v1).cross(v2 - v1);
    Eigen::Vector3d x30 = q3start - q1start;
    Eigen::Vector3d v30 = v3 - v1;
    planePoly3D(x10, x20, x30, v10, v20, v30, e1);

    if(e1.empty())
        return false;

    x10 = q0start - q2start;
    v10 = v0 - v2;
    x20 = (q1start - q2start).cross(q3start - q2start);
    v20 = (v1 - v2).cross(v3 - v2);
    x30 = q1start - q2start;
    v30 = v1 - v2;
    planePoly3D(x10, x20, x30, v10, v20, v30, e2);

    if(e2.empty())
        return false;

    x10 = q0start - q3start;
    v10 = v0 - v3;
    x20 = (q2start - q3start).cross(q1start - q3start);
    v20 = (v2 - v3).cross(v1 - v3);
    x30 = q2start - q3start;
    v30 = v2 - v3;
    planePoly3D(x10, x20, x30, v10, v20, v30, e3);

    if(e3.empty())
        return false;

    x10 = q0start - q1start;
    x20 = q2start - q1start;
    x30 = q3start - q1start;
    v10 = v0 - v1;
    v20 = v2 - v1;
    v30 = v3 - v1;
    distancePoly3D(x10, x20, x30, v10, v20, v30, mind, coplane);

    if(coplane.empty())
        return false;

    bool col = false;
    double mint = 1.0;
    for (int i = 0; i < (int) coplane.size(); i++)
    {
        for (int j = 0; j < (int) e1.size(); j++)
        {
            for (int k = 0; k < (int) e2.size(); k++)
            {
                for (int l = 0; l < (int) e3.size(); l++)
                {
                    vector<TimeInterval> intervals;
                    intervals.push_back(coplane[i]);
                    intervals.push_back(e1[j]);
                    intervals.push_back(e2[k]);
                    intervals.push_back(e3[l]);
                    if(TimeInterval::overlap(intervals))
                    {
                        mint = std::min(TimeInterval::intersect(intervals).l, mint);
                        col = true;
                    }
                }
            }
        }
    }

    if(col)
    {
        t = mint;
        return true;
    }
    return false;
}
