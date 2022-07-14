//
// Created by aldof on 4/10/2022.
//

#include <memory>
#include "CXPBDContinuousCollisionDetection.h"
#include "../../shared/rpoly.h"
#include "CXPBDAABB.h"

using namespace std;

static void Barycentric(float out[3], const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& Q)
{
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d QA = A - Q;
    Eigen::Vector3d QB = B - Q;

    float divisor = AB.dot(AB);

    out[0] = QB.dot(AB);
    out[1] = -QA.dot(AB);
    out[2] = divisor;
}

// Convert a point Q from Cartesian coordinates to Barycentric coordinates (u, v, w)
// with respect to a triangle ABC.
// The last output value is the divisor.
static void Barycentric(float out[4], const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C,
                        const Eigen::Vector3d& Q)
{
    Eigen::Vector3d AB = B - A;
    Eigen::Vector3d AC = C - A;

    Eigen::Vector3d QA = A - Q;
    Eigen::Vector3d QB = B - Q;
    Eigen::Vector3d QC = C - Q;

    Eigen::Vector3d QB_x_QC = QB.cross(QC);
    Eigen::Vector3d QC_x_QA = QC.cross(QA);
    Eigen::Vector3d QA_x_QB = QA.cross(QB);

    Eigen::Vector3d AB_x_AC = AB.cross(AC);

    float divisor = AB_x_AC.dot(AB_x_AC);

    out[0] = QB_x_QC.dot(AB_x_AC);
    out[1] = QC_x_QA.dot(AB_x_AC);
    out[2] = QA_x_QB.dot(AB_x_AC);
    out[3] = divisor;
}


// Return the closest point on a triangle ABC to a point Q.
bool SolvePoint(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C, const Eigen::Vector3d& Q,
                Eigen::Vector3d& p)
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
    TimeInterval isect(0.0, std::numeric_limits<double>::infinity());
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
    t1 = min(std::numeric_limits<double>::infinity(), t1);
    t2 = min(std::numeric_limits<double>::infinity(), t2);

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
            intervals.push_back(TimeInterval(0, std::numeric_limits<double>::infinity()));
        return;
    }

    // check intervals
    if (roots > 0)
    {
        std::sort(time, time + roots);
        if (time[0] >= 0)
            checkInterval(0, time[0], op, reducedDegree, intervals, pos);
        for (int i = 0; i < roots - 1; i++) {
            // if (!((time[i] < 0 && time[i + 1] < 0) || (time[i] > 1.0 && time[i + 1] > 1.0)))
            if (!((time[i] < 0 && time[i + 1] < 0) ))
                checkInterval(time[i], time[i + 1], op, reducedDegree, intervals, pos);
        }
        if (time[roots - 1] <= 1.0)
            checkInterval(time[roots - 1], std::numeric_limits<double>::infinity(), op, reducedDegree, intervals, pos);
    }
    else
    {
        checkInterval(0.0, std::numeric_limits<double>::infinity(), op, reducedDegree, intervals, pos);
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

// Dynamic collison : Checks for collison along path of goal and proxy
bool findCollisions(Eigen::Vector3d& goal, Eigen::Vector3d& proxy, double toolRadius,
                       cXPBDDeformableMesh* model, std::vector<ColInfo*>& collisions)
{

    // flag
    bool ret = 0;

    // Builds a bounding box for the tool
    auto toolBB_ = buildAABB(goal,proxy,toolRadius);

    // Get the object primitives
    const auto& f_ = model->faces();
    const auto& v_ = model->positions();
    const auto& vlast_ = model->positions();
    auto& N_ = model->normals();
    auto& Nlast_ = model->normals_last();

    std::vector<Collision> potentialCollisions;
    intersect(toolBB_, model->bb(), potentialCollisions);
    int it_count = 0;

    for (const auto& it : potentialCollisions) {
        // define collision info
        Eigen::Vector3i face = f_.row(it.collidingTriangle2);
        Eigen::Vector3d A0 = vlast_.row(face(0));
        Eigen::Vector3d B0 = vlast_.row(face(1));
        Eigen::Vector3d C0  = vlast_.row(face(2));
        Eigen::Vector3d N0 = Nlast_.row(it.collidingTriangle2);
        Eigen::Vector3d A1 = v_.row(face(0));
        Eigen::Vector3d B1 = v_.row(face(1));
        Eigen::Vector3d C1 = v_.row(face(2));
        Eigen::Vector3d N1 = N_.row(it.collidingTriangle2);

        // time of collision
        double t_;

        // index counter
        int counter = 0;

        if (CTCD::vertexFaceCTCD(
                proxy ,
                A0.transpose(),
                B0.transpose(),
                C0.transpose(),
                goal ,
                A1.transpose(),
                B1.transpose(),
                C1.transpose(),
                1e-6,
                t_))
        {
            // check what type of collision it is

            // Create a new collision object
            collisions.emplace_back(new ColInfo);

            // update collision info
            collisions[counter]->t = t_;

            // get the triangle on collision
            Eigen::Vector3d Ac;
            Eigen::Vector3d Bc;
            Eigen::Vector3d Cc;
            Eigen::Vector3d Nc;

            // get the proxy on collision
            Eigen::Vector3d Pc;

            // save normals
            collisions[counter]->normal0 = N0;
            collisions[counter]->normal1 = N1;

            // case where the proxy collides within the time interval
            if (t_ < 1.0)
            {

                // triangle vertices and normal on collision
                Ac = (A1 - A0) * t_ + A0;
                Bc = (B1 - B0) * t_ + B0;
                Cc = (C1 - C0) * t_ + C0;
                Nc = (N1 - N0) * t_ + N0;

                // proxy position on collision
                Pc = (goal - proxy) * t_ + proxy;

                // get barycentric coordinates
                double area = ((Ac - Bc).cross(Ac - Cc)).norm()/2;
                double alpha = ((Pc - Bc).cross(Pc - Cc)).norm()/(area*2);
                double beta = ((Pc - Cc).cross(Pc - Ac)).norm()/(area*2);
                double gamma = 1 - alpha - beta;

                // std::cout << " : alpha , " << alpha << " : beta , " << beta << " : gamma , " << gamma << std::endl;

                // check if center of proxy is within triangle
                if ( 0 <= alpha <= 1 && 0 <= beta <= 1 && 0 <= gamma <= 1)
                {
                    // if inside then face collision
                    collisions[counter]->type = FACECOLLISION;
                    collisions[counter]->alpha = alpha;
                    collisions[counter]->beta = beta;
                    collisions[counter]->gamma = gamma;
                    collisions[counter]->normalc = Nc;

                    // std::cout << "INSIDE" << std::endl;

                }
                else
                {
                    // if not inside then find the closest point on triangle
                    // and test is closer than tool radius
                    // std::cout << "OUTSIDE" << std::endl;
                }

            }
            // proxy collides outside of time interval
            else
            {
                // get barycentric coordinates
                double area = ((A1 - B1).cross(A1 - C1)).norm()/2;
                double alpha = ((proxy - B1).cross(proxy - C1)).norm()/(area*2);
                double beta = ((proxy - C1).cross(proxy - A1)).norm()/(area*2);
                double gamma = 1 - alpha - beta;

                // check if center of proxy is within triangle
                if ( 0 <= alpha <= 1 && 0 <= beta <= 1 && 0 <= gamma <= 1)
                {
                    // if inside then face collision
                    collisions[counter]->type = FACECOLLISION;
                    collisions[counter]->alpha = alpha;
                    collisions[counter]->beta = beta;
                    collisions[counter]->gamma = gamma;
                    collisions[counter]->normalc = Nc;

                    // Find closest point on the triangle
                    Eigen::Vector3d cp;

                    // find the closest point on the triangle
                    SolvePoint(A1,B1,C1,proxy,cp);

                    // TODO: MAKE SURE TO CORRECTLY HANDLE COLLISION WITH SPHERE

                    if ((proxy - cp).squaredNorm() < toolRadius*toolRadius)
                    {

                    }

                    // std::cout << "INSIDE" << std::endl;

                }
                else
                {
                    // if not inside then find the closest point on triangle
                    // and test is closer than tool radius
                    // std::cout << "OUTSIDE" << std::endl;
                }

            }

            // flag
            ret = 1;

            // increment counter
            counter++;
        }

    }
    return ret;

}

// Static collison : Just checks for collision with proxy
bool findCollisions(Eigen::Vector3d& proxy, double toolRadius,
                    cXPBDDeformableMesh* model, std::vector<ColInfo*>& collisions)
{

    // flag
    bool ret = 0;

    // Builds a bounding box for the tool
    auto toolBB_ = buildAABB(proxy,toolRadius);

    // Get the object primitives
    const auto& f_ = model->faces();
    const auto& v_ = model->positions();
    const auto& vlast_ = model->positions();
    auto& N_ = model->normals();
    auto& Nlast_ = model->normals_last();

    std::vector<Collision> potentialCollisions;
    intersect(toolBB_, model->bb(), potentialCollisions);
    int it_count = 0;

    for (const auto& it : potentialCollisions) {
        // define collision info
        Eigen::Vector3i face = f_.row(it.collidingTriangle2);
        Eigen::Vector3d A0 = vlast_.row(face(0));
        Eigen::Vector3d B0 = vlast_.row(face(1));
        Eigen::Vector3d C0  = vlast_.row(face(2));
        Eigen::Vector3d N0 = Nlast_.row(it.collidingTriangle2);
        Eigen::Vector3d A1 = v_.row(face(0));
        Eigen::Vector3d B1 = v_.row(face(1));
        Eigen::Vector3d C1 = v_.row(face(2));
        Eigen::Vector3d N1 = N_.row(it.collidingTriangle2);

        // time of collision
        double t_;

        // index counter
        int counter = 0;

        if (CTCD::vertexFaceCTCD(
                proxy ,
                A0.transpose(),
                B0.transpose(),
                C0.transpose(),
                A1.transpose(),
                B1.transpose(),
                C1.transpose(),
                1e-6,
                t_))
        {
            // check what type of collision it is

            // Create a new collision object
            collisions.emplace_back(new ColInfo);

            // update collision info
            collisions[counter]->t = t_;

            // get the triangle on collision
            Eigen::Vector3d Ac;
            Eigen::Vector3d Bc;
            Eigen::Vector3d Cc;
            Eigen::Vector3d Nc;

            // get the proxy on collision
            Eigen::Vector3d Pc;

            // save normals
            collisions[counter]->normal0 = N0;
            collisions[counter]->normal1 = N1;

            std::cout << t_ << std::endl;

            // case where the proxy collides within the time interval
            if (t_ < 1.0)
            {

                // triangle vertices and normal on collision
                Ac = (A1 - A0) * t_ + A0;
                Bc = (B1 - B0) * t_ + B0;
                Cc = (C1 - C0) * t_ + C0;
                Nc = (N1 - N0) * t_ + N0;

                // proxy position on collision
                //Pc = (goal - proxy) * t_ + proxy;

                // get barycentric coordinates
                double area = ((Ac - Bc).cross(Ac - Cc)).norm()/2;
                double alpha = ((Pc - Bc).cross(Pc - Cc)).norm()/(area*2);
                double beta = ((Pc - Cc).cross(Pc - Ac)).norm()/(area*2);
                double gamma = 1 - alpha - beta;

                // std::cout << " : alpha , " << alpha << " : beta , " << beta << " : gamma , " << gamma << std::endl;

                // check if center of proxy is within triangle
                if ( 0 <= alpha <= 1 && 0 <= beta <= 1 && 0 <= gamma <= 1)
                {
                    // if inside then face collision
                    collisions[counter]->type = FACECOLLISION;
                    collisions[counter]->alpha = alpha;
                    collisions[counter]->beta = beta;
                    collisions[counter]->gamma = gamma;
                    collisions[counter]->normalc = Nc;

                    // std::cout << "INSIDE" << std::endl;

                }
                else
                {
                    // if not inside then find the closest point on triangle
                    // and test is closer than tool radius
                    // std::cout << "OUTSIDE" << std::endl;
                }

            }
                // proxy collides outside of time interval
            else
            {
                //std::cout << "COLLISION t > 1 " << std::endl;
                Eigen::Vector3d cp;

                SolvePoint(A1,B1,C1,proxy,cp);
                std::cout << cp.transpose() << std::endl;
            }

            // flag
            ret = 1;

            // increment counter
            counter++;
        }

    }
    return ret;

}

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
                        mint = TimeInterval::intersect(intervals).l;
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

bool CTCD::vertexFaceCTCD(const Eigen::Vector3d& q0start,
                          const Eigen::Vector3d& q1start,
                          const Eigen::Vector3d& q2start,
                          const Eigen::Vector3d& q3start,
                          const Eigen::Vector3d& q1end,
                          const Eigen::Vector3d& q2end,
                          const Eigen::Vector3d& q3end,
                          double eta,
                          double& t)
{

    double mind = eta*eta;

    Eigen::Vector3d v0 = q0start;
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
                        mint = TimeInterval::intersect(intervals).l;
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
