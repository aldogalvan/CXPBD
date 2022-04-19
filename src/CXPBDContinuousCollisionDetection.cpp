//
// Created by aldof on 4/10/2022.
//

#include <memory>
#include "CXPBDContinuousCollisionDetection.h"
#include "../shared/rpoly.h"

using namespace std;

ColInfo findCollisions(cXPBDDeformableMesh* def_object,
                       cXPBDTool* tool)
{
    ColInfo ret;
    auto& vertCollisions = ret.vertCollisions;
    auto& faceCollisions = ret.faceCollisions;

    vertCollisions.resize(tool->numVerts());
    faceCollisions.resize(tool->numFaces());

    def_object->buildAABBBoundaryBox();
    tool->buildAABBBoundaryBox();

    const auto& F_def = def_object->faces();
    const auto& F_tool = tool->faces();
    const auto& V_def =def_object->positions();
    const auto& V_tool = tool->positions();
    const auto& V_def_last = def_object->positions_last();
    const auto& V_tool_last = tool->positions_last();
    std::vector<Collision> potentialCollisions;

    intersect(tool->bb(), def_object->bb(), potentialCollisions);

    for (const auto& it : potentialCollisions) {
        Eigen::Vector3i face = F_def.row(it.collidingTriangle2);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = F_tool(it.collidingTriangle1, vert);
            double t;

            if (CTCD::vertexFaceCTCD(
                    V_tool_last.row(vidx).transpose(),
                    V_def_last.row(face(0)).transpose(),
                    V_def_last.row(face(1)).transpose(),
                    V_def_last.row(face(2)).transpose(),
                    V_tool.row(vidx).transpose(),
                    V_def.row(face(0)).transpose(),
                    V_def.row(face(1)).transpose(),
                    V_def.row(face(2)).transpose(),
                    1e-6,
                    t))
            {
                vertCollisions[vidx].emplace(it.collidingTriangle2);
            }
        }
    }

    for (const auto& it : potentialCollisions) {
        Eigen::Vector3i face = F_tool.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = F_def(it.collidingTriangle2, vert);
            double t;
            // TODO Implement Vertex-Face CTCD
            if (CTCD::vertexFaceCTCD(
                    V_def_last.row(vidx).transpose(),
                    V_tool_last.row(face(0)).transpose(),
                    V_tool_last.row(face(1)).transpose(),
                    V_tool_last.row(face(2)).transpose(),
                    V_def.row(vidx).transpose(),
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
    CTCD fun;

    Eigen::Vector3d ai = q1start - q0start;
    Eigen::Vector3d bi = q2start - q0start;
    Eigen::Vector3d ci = q3start - q0start;
    Eigen::Vector3d af = q1end - q0end;
    Eigen::Vector3d bf = q2end - q0end;
    Eigen::Vector3d cf = q3end - q0end;

    auto [t_0, t_1, t_2, t_3] = fun.GetTerms(ai, bi, ci, af, bf, cf);


    int size;
    int deg;

    if (t_3 <= 1e-6 && t_3 >= -1e-6) {
        size = 3;
        deg = 2;
        if (t_2 <= 1e-6 && t_2 >= -1e-6) {
            size = 2;
            deg = 1;
            if (t_1 <= 1e-6 && t_1 >= -1e-6){
                return false;
            }
        }
    }
    else{
        size = 4;
        deg = 3;
    }

    double poly[size];

    if (size == 4) {
        poly[3] = t_3;
        poly[2] = t_2;
        poly[1] = t_1;
        poly[0] = t_0;
    }
    else if (size == 3){
        poly[2] = t_2;
        poly[1] = t_1;
        poly[0] = t_0;
    }
    else{
        poly[1] = t_1;
        poly[0] = t_0;
    }

    RootFinder rf;
    int sol;

    double real[4] = {};
    double im[4] = {};

    sol = rf.rpoly( poly, deg, real, im);
    if (sol == -1)
    {
        return false;
    }
    for (int i = 0; i < sol; i++) {
        if ((real[i] >= 0) && (real[i] <= 1)) {
            auto pos = q0start + (q0end - q0start) * real[i];
            auto p1 = q1start + (q1end - q1start) * real[i];
            auto p2 = q2start + (q2end - q2start) * real[i];
            auto p3 = q3start + (q3end - q3start) * real[i];
            float u;
            float v;
            float w;
            fun.Barycentric(pos, p1, p2, p3, u, v, w);
            if (u >= 0 && u <= 1 && v >= 0 && v <= 1 && w >= 0 && w <= 1) {
                return true;

            }
        }
    }
    return false;
}


std::tuple<double, double, double, double>
CTCD::GetTerms( Eigen::Vector3d& x0, Eigen::Vector3d& y0, Eigen::Vector3d& z0, Eigen::Vector3d& x,
                Eigen::Vector3d& y, Eigen::Vector3d& z){
    double c = x0.cross(y0).dot(z0);
    double t = x0.cross(y0).dot(z - z0) + x0.cross(y - y0).dot(z0) + (x - x0).cross(y0).dot(z0);
    double t_2 = x0.cross(y - y0).dot(z - z0) + (x - x0).cross(y0).dot(z - z0) + (x - x0).cross(y - y0).dot(z0);
    double t_3 = (x - x0).cross(y - y0).dot(z - z0);

    return std::make_tuple(c, t, t_2, t_3);
}

void CTCD::Barycentric(Eigen::Vector3d p, Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, float &u, float &v, float &w)
{
    Eigen::Vector3d v0 = b - a, v1 = c - a, v2 = p - a;
    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;
    v = (d11 * d20 - d01 * d21) / denom;
    w = (d00 * d21 - d01 * d20) / denom;
    u = 1.0f - v - w;
}
