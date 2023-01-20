#include "collision.h"
#include "collision_kernel.cuh"
#include "../../shared/rpoly.h"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

/////////////////////////////////////////////////////////////////
/// MAIN FUNCTION TO FIND COLLISIONS
/////////////////////////////////////////////////////////////////


void findCollisions(const float* d_goal, float* d_proxy, meshObject* obj)
{

    // broad phase collision
    auto potCol = CTCD::broadPhase(d_goal, d_proxy, obj);
    // narrow phase collision detection
    //CTCD::narrowPhase(d_goal,d_proxy,obj,potCol);

}


static void Barycentric(double out[3], const Vector3d& A, const Vector3d& B, const Vector3d& Q)
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
static void Barycentric(double out[4], const Vector3d& A, const Vector3d& B, const Vector3d& C,
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

    double divisor = AB_x_AC.dot(AB_x_AC);

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
    double wAB[3], wBC[3], wCA[3];
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
    double wABC[4];
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
    TimeInterval isect(0.0, 1.0);
    for(std::vector<TimeInterval>::const_iterator it = intervals.begin(); it != intervals.end(); ++it)
    {
        isect.l = max(it->l, isect.l);
        isect.u = min(it->u, isect.u);
    }
    return isect;
}

vector<int> CTCD::broadPhase(const float* d_goal, float* d_proxy, meshObject* obj)
{

    /*
    // Builds a bounding box for the tool
    aabb* rhs = new aabb;

    // initializes the bounding box
    rhs->lower.x  = numeric_limits<float>::infinity();
    rhs->lower.y  = numeric_limits<float>::infinity();
    rhs->lower.z  = numeric_limits<float>::infinity();

    rhs->upper.x  = -numeric_limits<float>::infinity();
    rhs->upper.y  = -numeric_limits<float>::infinity();
    rhs->upper.z  = -numeric_limits<float>::infinity();

    // build the bounding box
    rhs->lower.x = min(rhs->lower.x, d_goal[0]);
    rhs->upper.x = max(rhs->upper.x, d_goal[0]);
    rhs->lower.x = min(rhs->lower.x, d_proxy[0]);
    rhs->upper.x = max(rhs->upper.x, d_proxy[0]);

    rhs->lower.y = min(rhs->lower.y, d_goal[1]);
    rhs->upper.y = max(rhs->upper.y, d_goal[1]);
    rhs->lower.y = min(rhs->lower.y, d_proxy[1]);
    rhs->upper.y = max(rhs->upper.y, d_proxy[1]);

    rhs->lower.z = min(rhs->lower.z, d_goal[2]);
    rhs->upper.z = max(rhs->upper.z, d_goal[2]);
    rhs->lower.z = min(rhs->lower.z, d_proxy[2]);
    rhs->upper.z = max(rhs->upper.z, d_proxy[2]);

    // the device pointer
    aabb* drhs;
    cudaMalloc((void**)&drhs,2*sizeof(float3));
    cudaMemcpy(drhs,rhs,2*sizeof(float3),cudaMemcpyHostToDevice);

    // the vector containing the possible collisions
    int* dcol;
    cudaMalloc((void**)&dcol,obj->nfaces*sizeof(int));

    // TODO: CHECK IF ASYNCHRONOUS COLLISION DETECTION IS BETTER
    // BASICALLY ONLY COMPUTE THE BOUNDING BOX WHEN DYNAMICS ARE DONE
    // build bounding box
    aabb** dlhs;
    cudaMalloc((void**)&dlhs,2*obj->nfaces*sizeof(float3));
    simple_bounding_box<<<2,1024>>>(dlhs,obj->d_fi,obj->d_xlast,obj->d_x,obj->nfaces);

    // computes the lhs
    simple_broad_phase<<<2,1024>>>(dlhs,drhs,dcol,obj->nfaces);

    // return
    //thrust::device_vector<int> d_pcol(dcol, dcol + obj->nfaces);
    //thrust::host_vector<int> h_pcol = d_pcol;
    std::vector<int> ret;


    // delete initialized pointer
    delete rhs;
    cudaFree(drhs);
    cudaFree(dlhs);
    cudaFree(dcol);
     */

    return ret;
}


void CTCD::narrowPhase(const float* d_goal, const float* d_proxy, const meshObject* obj,
                       vector<int>& potentialCollisions)
{

    /*
    int ncol = potentialCollisions.size();
    int* collisions = potentialCollisions.data();
    int* d_collisions;
    float* d_goal; float* d_proxy;
    float* t_;

    cudaMalloc((void**)collisions,ncol*sizeof(int));
    cudaMalloc((void**)t_,ncol*sizeof(float));
    cudaMalloc((void**)d_goal,3*sizeof(float));
    cudaMalloc((void**)d_proxy,3*sizeof(float));

    cudaMemcpy(d_collisions,collisions,ncol*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(d_goal,goal,3*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(d_proxy,proxy,3*sizeof(float),cudaMemcpyHostToDevice);
    */
    /*
    for (const auto& it : potentialCollisions) {

        // define collision info
        Vector3d face = f_.row(it);
        Vector3d A0 = plast_.row(face(0));
        Vector3d B0 = plast_.row(face(1));
        Vector3d C0 = plast_.row(face(2));
        Vector3d N0 = -Nlast_.row(it.collidingTriangle2).normalized();
        Vector3d A1 = p_.row(face(0));
        Vector3d B1 = p_.row(face(1));
        Vector3d C1 = p_.row(face(2));
        Vector3d N1 = -N_.row(it.collidingTriangle2).normalized();

        // time of collision
        double t_;

        bool flag = 0;

        if (CTCD::vertexFaceCTCD(proxy,A0.transpose(),B0.transpose(),C0.transpose(),
                                 goal,A1.transpose(),B1.transpose(),C1.transpose(),1e-6,t_))
        {

            double ep = 1e-6;

            //! TODO: Separate the collision handling from collision detection

            // Create a new collision object
            ColInfo *col = new ColInfo;

            // update collision info
            col->t = t_;
            col->triangle = face;

            // Triangle at collision
            Vector3d Ac = A0 + t_*(A1 - A0);
            Vector3d Bc = B0 + t_*(B1 - B0);
            Vector3d Cc = C0 + t_*(C1 - C0);

            // Normal at collision
            Vector3d Nc = N0 + t_*(N1 - N0);

            // proxy at collision
            Vector3d proxy_c = proxy + t_*(goal - proxy);

            // add a threshold

            // Barycentric coordinates on collision
            double out[4];
            Barycentric(out,Ac,Bc,Cc,proxy_c);

            // move to barycentric coordinates
            proxy = out[0]*A1/out[3] + out[1]*B1/out[3] + out[2]*C1/out[3];

            // get the equation of the plane defined by the triangle
            double A, B, C, D;
            A = N1(0); B = N1(1); C = N1(2); D = -N1.dot(A1);

            // move the proxy tangent to the plane
            double dist = abs(A*goal(0) + B*goal(1) + C*goal(2) + D) / sqrt(pow(A,2) + pow(B,2) + pow(C,2));
            Vector3d dgoal = goal + dist*N1;
            proxy = dgoal;

            // check if the proxy and goal point are
            // on opposite sides of the triangle
            // if so then add collision
            if ( (proxy - A1).dot(N1) > 0 && (goal - A1).dot(N1) <= 0)
            {
                collisions.emplace_back(col);
            }
                // check if the proxy and goal point are
                // both on outside of triangle
                // if so then continue
            else if ( (proxy - A1).dot(N1) >= 0 && (goal - A1).dot(N1) >= 0)
            {
                continue;
            }
                // check if the proxy and goal point are
                // on opposite sides of triangle in degenerate case
                // if so then project the proxy to the goal and continue
            else if ((proxy - A1).dot(N1) < 0 && (goal - A1).dot(N1) >= 0)
            {
                proxy = goal;
                continue;
            }
                // check if the proxy and goal point are
                // both on inside of triangle
                // if so then project to surface and add
                // a error threshold
            else if((proxy - A1).dot(N1) <= 0 && (goal - A1).dot(N1) <= 0)
            {
                //! NEED TO PERFORM REPROJECTION
                double dist = abs(A*proxy(0) + B*proxy(1) + C*proxy(2) + D) / sqrt(pow(A,2) + pow(B,2) + pow(C,2));
                proxy -= (dist + ep)*N1;
                //proxy += 2*dist*N1;
                //proxy += N1 * ep;
                collisions.emplace_back(col);
                //cout << (proxy - A1).dot(N1) << endl;
                //cout << (proxy - A1).dot(N1) << endl;
                //cout << out[0]/out[3] << "," << out[1]/out[3] << "," << out[2]/out[3] << endl;
            }

            // Now we slide along the point subject to friction
            // very basic (will need to add)


            // adjust a bit in normal direction to account
            // for numerical errors
            /*
            if ((proxy - A1).dot(N1) < 1)
            {
                proxy += N1.normalized() * ep;
            }
            else
            {
                proxy -= N1.normalized() * ep;
            }
             */

/*
            ret = 1;
        }


        else if (CTCD::vertexEdgeCTCD(proxy,A0,B0,goal,A1,B1,1e-6,t_))
        {
            cout << "EDGE COLLISION 1" << endl;
        }
        else if (CTCD::vertexEdgeCTCD(proxy,B0,C0,goal,B1,C1,1e-6,t_))
        {
            cout << "EDGE COLLISION 2" << endl;
        }
        else if (CTCD::vertexEdgeCTCD(proxy,A0,C0,goal,A1,C1,1e-6,t_))
        {
            cout << "EDGE COLLISION 3" << endl;
        }
        else if (CTCD::vertexVertexCTCD(proxy ,A0,goal,A1 ,1E-6 ,t_))
        {
            cout << "VERTEX COLLISION 1" << endl;
        }
        else if (CTCD::vertexVertexCTCD(proxy,B0,goal,B1 ,1e-6 ,t_))
        {
            cout << "VERTEX COLLISION 2" << endl;
        }
        else if (CTCD::vertexVertexCTCD(proxy,C0,goal,C1 ,1e-6 ,t_))
        {
            cout << "VERTEX COLLISION 3" << endl;
        }
    }
    */


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

void CTCD::barycentricPoly3D(const Vector3d &x10,
                             const Vector3d &x20,
                             const Vector3d &x30,
                             const Vector3d &v10,
                             const Vector3d &v20,
                             const Vector3d &v30,
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

void CTCD::planePoly3D(const Vector3d &x10,
                       const Vector3d &x20,
                       const Vector3d &x30,
                       const Vector3d &v10,
                       const Vector3d &v20,
                       const Vector3d &v30,
                       vector<TimeInterval> &result)
{

    double op[4];
    op[0] = v10.dot(v20.cross(v30));
    op[1] = x10.dot(v20.cross(v30)) + v10.dot(x20.cross(v30)) + v10.dot(v20.cross(x30));
    op[2] = x10.dot(x20.cross(v30)) + x10.dot(v20.cross(x30)) + v10.dot(x20.cross(x30));
    op[3] = x10.dot(x20.cross(x30));
    findIntervals(op, 3, result, true);
}

void CTCD::distancePoly3D(const Vector3d &x10,
                          const Vector3d &x20,
                          const Vector3d &x30,
                          const Vector3d &v10,
                          const Vector3d &v20,
                          const Vector3d &v30,
                          double minDSquared,
                          vector<TimeInterval> &result)
{
    double A = v10.dot(v20.cross(v30));
    double B = x10.dot(v20.cross(v30)) + v10.dot(x20.cross(v30)) + v10.dot(v20.cross(x30));
    double C = x10.dot(x20.cross(v30)) + x10.dot(v20.cross(x30)) + v10.dot(x20.cross(x30));
    double D = x10.dot(x20.cross(x30));
    Vector3d E = x20.cross(x30);
    Vector3d F = x20.cross(v30) + v20.cross(x30);
    Vector3d G = v20.cross(v30);

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


bool CTCD::edgeEdgeCTCD(const Vector3d &q0start,
                        const Vector3d &p0start,
                        const Vector3d &q1start,
                        const Vector3d &p1start,
                        const Vector3d &q0end,
                        const Vector3d &p0end,
                        const Vector3d &q1end,
                        const Vector3d &p1end,
                        double eta,
                        double &t)
{
    double minD = eta * eta;

    // time intervals during which v is colinear with the edge, on the side of e1 towards e2, and on the side of e2 towards e1
    std::vector<TimeInterval> rawcoplane, a0, a1, b0, b1;

    Vector3d x10 = p0start - p1start;
    Vector3d x20 = p0start - q0start;
    Vector3d x30 = p1start - q1start;

    Vector3d vp0 = p0end-p0start;
    Vector3d vp1 = p1end-p1start;
    Vector3d vq0 = q0end-q0start;
    Vector3d vq1 = q1end-q1start;

    Vector3d v10 = vp0 - vp1;
    Vector3d v20 = vp0 - vq0;
    Vector3d v30 = vp1 - vq1;

    distancePoly3D(x10, x20, x30, v10, v20, v30, minD, rawcoplane);

    // check for parallel edges
    std::vector<TimeInterval> coplane;
    std::vector<TimeInterval> parallel;

    for(size_t i=0; i<rawcoplane.size(); i++)
    {
        double midt = (rawcoplane[i].u + rawcoplane[i].l) / 2;
        x10 = (q0start - p0start) + midt * (vq0 - vp0);
        x20 = (q1start - p1start) + midt * (vq1 - vp1);

        if (x10.cross(x20).norm() < 1e-8)
        {
            parallel.push_back(rawcoplane[i]);
        }
        else
        {
            coplane.push_back(rawcoplane[i]);
        }
    }

    if(coplane.empty())
        return false;

    x10 = p1start - q1start;
    v10 = vp1 - vq1;
    x20 = p0start - q0start;
    v20 = vp0 - vq0;
    x30 = q0start - q1start;
    v30 = vq0 - vq1;
    barycentricPoly3D(x10, x20, x30, v10, v20, v30, a0);
    if(a0.empty())
        return false;

    x20 = q0start - p0start;
    v20 = vq0 - vp0;
    x30 = p0start - q1start;
    v30 = vp0 - vq1;
    barycentricPoly3D(x10, x20, x30, v10, v20, v30, a1);
    if(a1.empty())
        return false;

    x10 = p0start - q0start;
    v10 = vp0 - vq0;
    x20 = p1start - q1start;
    v20 = vp1 - vq1;
    x30 = q1start - q0start;
    v30 = vq1 - vq0;
    barycentricPoly3D(x10, x20, x30, v10, v20, v30, b0);
    if(b0.empty())
        return false;

    //x10 = p0 - q0;
    //v10 = vp0 - vq0;
    x20 = q1start - p1start;
    v20 = vq1 - vp1;
    x30 = p1start - q0start;
    v30 = vp1 - vq0;

    barycentricPoly3D(x10, x20, x30, v10, v20, v30, b1);
    if(b1.empty())
        return false;

    // check intervals for overlap
    bool col = false;
    double mint = 1.0;
    for (int i = 0; i < (int) coplane.size(); i++)
    {
        for (int j = 0; j < (int) a0.size(); j++)
        {
            for (int k = 0; k < (int) a1.size(); k++)
            {
                for (int l = 0; l < (int) b0.size(); l++)
                {
                    for (int m = 0; m < (int) b1.size(); m++)
                    {
                        vector<TimeInterval> intervals;
                        intervals.push_back(coplane[i]);
                        intervals.push_back(a0[j]);
                        intervals.push_back(a1[k]);
                        intervals.push_back(b0[l]);
                        intervals.push_back(b1[m]);
                        if (TimeInterval::overlap(intervals) )
                        {
                            TimeInterval isect = TimeInterval::intersect(intervals);
                            bool skip = false;
                            for(int p = 0; p < (int)parallel.size(); p++)
                            {
                                vector<TimeInterval> pcheck;
                                pcheck.push_back(isect);
                                pcheck.push_back(parallel[p]);
                                if(TimeInterval::overlap(pcheck) )
                                {
                                    skip = true;
                                    break;
                                }
                            }
                            if(!skip)
                            {
                                mint = min(mint, isect.l);
                                col = true;
                            }
                        }
                    }
                }
            }
        }
    }

    // handle parallel edges
    for(int i=0; i < (int)parallel.size(); i++)
    {

    }


    if(col)
    {
        t = mint;
        return true;
    }
    return false;
}

bool CTCD::vertexFaceCTCD(const Vector3d& q0start,
                          const Vector3d& q1start,
                          const Vector3d& q2start,
                          const Vector3d& q3start,
                          const Vector3d& q0end,
                          const Vector3d& q1end,
                          const Vector3d& q2end,
                          const Vector3d& q3end,
                          double eta,
                          double& t)
{
    double mind = eta*eta;

    Vector3d v0 = q0end - q0start;
    Vector3d v1 = q1end - q1start;
    Vector3d v2 = q2end - q2start;
    Vector3d v3 = q3end - q3start;

    // time intervals during which v is colinear with the edge,
    // on the side of e1 towards e2, and on the side of e2 towards e1
    vector<TimeInterval> coplane, e1, e2, e3;

    Vector3d x10 = q0start - q1start;
    Vector3d v10 = v0 - v1;
    Vector3d x20 = (q3start - q1start).cross(q2start - q1start);
    Vector3d v20 = (v3 - v1).cross(v2 - v1);
    Vector3d x30 = q3start - q1start;
    Vector3d v30 = v3 - v1;
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

bool CTCD::vertexFaceCTCD(const Vector3d& q0start,
                          const Vector3d& q1start,
                          const Vector3d& q2start,
                          const Vector3d& q3start,
                          const Vector3d& q1end,
                          const Vector3d& q2end,
                          const Vector3d& q3end,
                          double eta,
                          double& t)
{

    double mind = eta*eta;

    Vector3d v0 = q0start;
    Vector3d v1 = q1end - q1start;
    Vector3d v2 = q2end - q2start;
    Vector3d v3 = q3end - q3start;

    // time intervals during which v is colinear with the edge,
    // on the side of e1 towards e2, and on the side of e2 towards e1
    vector<TimeInterval> coplane, e1, e2, e3;

    Vector3d x10 = q0start - q1start;
    Vector3d v10 = v0 - v1;
    Vector3d x20 = (q3start - q1start).cross(q2start - q1start);
    Vector3d v20 = (v3 - v1).cross(v2 - v1);
    Vector3d x30 = q3start - q1start;
    Vector3d v30 = v3 - v1;
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

bool CTCD::vertexEdgeCTCD(const Vector3d &q0start,
                          const Vector3d &q1start,
                          const Vector3d &q2start,
                          const Vector3d &q0end,
                          const Vector3d &q1end,
                          const Vector3d &q2end,
                          double eta,
                          double &t)
{
    double op[5];
    double minD = eta*eta;
    Vector3d v0 = q0end-q0start;
    Vector3d v1 = q1end-q1start;
    Vector3d v2 = q2end-q2start;

    // time intervals during which v is colinear with the edge, on the side of e1 towards e2, and on the side of e2 towards e1
    vector<TimeInterval> colin, e1, e2;

    Vector3d ab = q2start - q1start;
    Vector3d ac = q0start - q1start;
    Vector3d cb = q2start - q0start;
    Vector3d vab = v2 - v1;
    Vector3d vac = v0 - v1;
    Vector3d vcb = v2 - v0;

    double c = ab.dot(ac);
    double b = ac.dot(vab) + ab.dot(vac);
    double a = vab.dot(vac);
    op[0] = a;
    op[1] = b;
    op[2] = c;
    findIntervals(op, 2, e1, true);
    if(e1.empty())
        return false;

    c = ab.dot(cb);
    b = cb.dot(vab) + ab.dot(vcb);
    a = vab.dot(vcb);

    op[0] = a;
    op[1] = b;
    op[2] = c;
    findIntervals(op, 2, e2, true);
    if(e2.empty())
        return false;

    double A = ab.dot(ab);
    double B = 2 * ab.dot(vab);
    double C = vab.dot(vab);
    double D = ac.dot(ac);
    double E = 2 * ac.dot(vac);
    double F = vac.dot(vac);
    double G = ac.dot(ab);
    double H = vab.dot(ac) + vac.dot(ab);
    double I = vab.dot(vac);
    op[4] = A * D - G * G - minD * A;
    op[3] = B * D + A * E - 2 * G * H - minD * B;
    op[2] = B * E + A * F + C * D - H * H - 2 * G * I - minD * C;
    op[1] = B * F + C * E - 2 * H * I;
    op[0] = C * F - I * I;
    findIntervals(op, 4, colin, false);
    if(colin.empty())
        return false;

    double mint = 1.0;
    bool col = false;
    for (int i = 0; i < (int) colin.size(); i++)
    {
        for (int j = 0; j < (int) e1.size(); j++)
        {
            for (int k = 0; k < (int) e2.size(); k++)
            {
                vector<TimeInterval> intervals;
                intervals.push_back(colin[i]);
                intervals.push_back(e1[j]);
                intervals.push_back(e2[k]);
                if(TimeInterval::overlap(intervals))
                {
                    mint = std::min(TimeInterval::intersect(intervals).l, mint);
                    col = true;
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

bool CTCD::vertexVertexCTCD(const Vector3d &q1start,
                            const Vector3d &q2start,
                            const Vector3d &q1end,
                            const Vector3d &q2end,
                            double eta, double &t)
{
    int roots = 0;
    double min_d = eta*eta;
    double t1 = 0, t2 = 0;
    Vector3d v1 = q1end-q1start;
    Vector3d v2 = q2end-q2start;

    // t^2 term
    double a = v1.dot(v1) + v2.dot(v2) - 2 * v1.dot(v2);
    // t term
    double b = 2 * (v1.dot(q1start) - v2.dot(q1start) - v1.dot(q2start) + v2.dot(q2start));
    // current distance - min_d
    double c = q1start.dot(q1start) + q2start.dot(q2start) - 2 * q1start.dot(q2start) - min_d;
    if (a != 0)
    {
        roots = getQuadRoots(a, b, c, t1, t2);
    }
    else if (b != 0)
    {
        t1 = -c / b;
        roots = 1;
    }
    else
    {
        if(c<=0)
        {
            t = 0;
            return true;
        }
        return false;
    }

    double op[3];
    op[0] = a;
    op[1] = b;
    op[2] = c;
    vector<TimeInterval> interval;

    if (roots == 2)
    {
        checkInterval(0, t1, op, 2, interval, false);
        if(!interval.empty())
        {
            t = 0;
            return true;
        }
        checkInterval(t1, t2, op, 2, interval, false);
        if(!interval.empty())
        {
            t = t1;
            return true;
        }
        checkInterval(t2, 1.0, op, 2, interval, false);
        if(!interval.empty())
        {
            t = t2;
            return true;
        }
        return false;
    }
    else if (roots == 1)
    {
        checkInterval(0, t1, op, 2, interval, false);
        if(!interval.empty())
        {
            t = 0;
            return true;
        }
        checkInterval(t1, 1.0, op, 2, interval, false);
        if(!interval.empty())
        {
            t = t1;
            return true;
        }
        return false;
    }
    checkInterval(0, 1.0, op, 2, interval, false);
    if(!interval.empty())
    {
        t = 0;
        return true;
    }
    return false;
}