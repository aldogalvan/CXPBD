//
// Created by aldo on 6/20/22.
//


#include "constraints/CXPBDEdgeLengthConstraint.h"
#include "constraints/CXPBDVolumeConstraint.h"
#include "constraints/CXPBDNeoHookeanConstraint.h"
#include "constraints/CXPBDBendingConstraint.h"
#include "constraints/CXPBDIsometricBendingConstraint.h"
#include "include/igl/per_face_normals.h"
#include "CXPBDTriangleMesh.h"

void cXPBDTriangleMesh::connectToChai3d(void)
{

    for (int i = 0; i < num_vertices; i++)
    {
        Eigen::Vector3d vert = p_.row(i);
        newVertex(vert.x(),vert.y(),vert.z());
    }

    for (int i = 0; i < num_faces; i++)
    {
        Eigen::Vector3i triangle = F_.row(i);
        newTriangle(triangle.x(),triangle.y(),triangle.z());
    }
}

void cXPBDTriangleMesh::updateChai3d(Eigen::MatrixXd& a_pos)
{
    for (int i = 0; i < num_vertices; i++)
    {
        Eigen::Vector3d vert = a_pos.row(i);
        m_vertices->setLocalPos(i,vert.x(),vert.y(),vert.z());
    }
}

Eigen::Vector3d cXPBDTriangleMesh::computeCentroid(void)
{
    Eigen::Vector3d sum;
    sum.setZero();

    for (int i = 0 ; i < num_vertices ; i++)
    {
        sum += p0_.row(i);
    }

    sum /= num_vertices;

    return sum;
}

void cXPBDTriangleMesh::setLocalPos(Eigen::Vector3d a_pos)
{
    Eigen::Vector3d centroid = computeCentroid();
    Eigen::Vector3d trans = a_pos - centroid;

    for (int i = 0 ; i < num_vertices ; i++)
    {
        p0_.row(i) += trans;
        p_.row(i) += trans;
    }

}

void cXPBDTriangleMesh::scaleObject(double a_scale)
{
    p0_ *= a_scale;
    plast_ *= a_scale;
    p_ *= a_scale;
}

void cXPBDTriangleMesh::constrain_edge_lengths(scalar_type const compliance , scalar_type const damping)
{
    auto const& positions = this->p0();
    auto const& E  = this->edges();


    for (auto i = 0u; i < E.rows(); ++i)
    {
        auto const edge = E.row(i);
        auto const e0   = edge(0);
        auto const e1   = edge(1);

        auto constraint = std::make_unique<cXPBDEdgeConstraint>(
                std::initializer_list<std::uint32_t>{
                        static_cast<std::uint32_t>(e0),
                        static_cast<std::uint32_t>(e1)},
                positions,
                compliance,
                damping);

        this->constraints().push_back(std::move(constraint));
    }
}

void cXPBDTriangleMesh::constrain_hinge_bending(const scalar_type compliance , scalar_type const damping)
{
    auto const& positions = this->p0();
    auto const& F = this->faces();
    std::map<std::pair<int, int>, std::vector<int>> edgemap;
    int nfaces = F.rows();

    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int nextj = (j + 1) % 3;
            int v1 = F(i, j);
            int v2 = F(i, nextj);
            if (v1 > v2)
                std::swap(v1, v2);
            edgemap[std::pair<int, int>(v1, v2)].push_back(i);
        }
    }

    int nhinges = 0;
    for (auto it : edgemap) {
        if (it.second.size() == 2)
            nhinges++;
    }
    H_.resize(nhinges, 4);
    int idx = 0;
    for (auto it : edgemap) {
        if (it.second.size() != 2)
            continue;
        std::set<int> hingeverts;
        for (int j = 0; j < 3; j++) {
            hingeverts.insert(F(it.second[0], j));
            hingeverts.insert(F(it.second[1], j));
        }
        int colidx = 0;
        for (auto v : hingeverts) {
            H_(idx, colidx) = v;
            colidx++;
        }
        idx++;
    }

    auto elements = H_;
    for (auto i = 0u; i < elements.rows(); ++i)
    {
        auto const element = elements.row(i);
        auto constraint    = std::make_unique<cXPBDBendingConstraint>(
                std::initializer_list<std::uint32_t>{
                        static_cast<std::uint32_t>(element(0)),
                        static_cast<std::uint32_t>(element(1)),
                        static_cast<std::uint32_t>(element(2)),
                        static_cast<std::uint32_t>(element(3))},
                positions,
                compliance,
                damping);

        this->constraints().push_back(std::move(constraint));
    }
}

void cXPBDTriangleMesh::constrain_isometric_bending(const scalar_type compliance, const scalar_type damping)
{
    auto const& positions = this->p0();
    auto const& F = this->faces();
    std::map<std::pair<int, int>, std::vector<int>> edgemap;
    int nfaces = F.rows();

    for (int i = 0; i < nfaces; i++) {
        for (int j = 0; j < 3; j++) {
            int nextj = (j + 1) % 3;
            int v1 = F(i, j);
            int v2 = F(i, nextj);
            if (v1 > v2)
                std::swap(v1, v2);
            edgemap[std::pair<int, int>(v1, v2)].push_back(i);
        }
    }

    int nhinges = 0;
    for (auto it : edgemap) {
        if (it.second.size() == 2)
            nhinges++;
    }
    H_.resize(nhinges, 4);
    int idx = 0;
    for (auto it : edgemap) {
        if (it.second.size() != 2)
            continue;
        std::set<int> hingeverts;
        for (int j = 0; j < 3; j++) {
            hingeverts.insert(F(it.second[0], j));
            hingeverts.insert(F(it.second[1], j));
        }
        int colidx = 0;
        for (auto v : hingeverts) {
            H_(idx, colidx) = v;
            colidx++;
        }
        idx++;
    }

    auto elements = H_;
    for (auto i = 0u; i < elements.rows(); ++i)
    {
        auto const element = elements.row(i);
        auto constraint    = std::make_unique<cXPBDIsometricBendingConstraint>(
                std::initializer_list<std::uint32_t>{
                        static_cast<std::uint32_t>(element(0)),
                        static_cast<std::uint32_t>(element(1)),
                        static_cast<std::uint32_t>(element(2)),
                        static_cast<std::uint32_t>(element(3))},
                positions,
                compliance,
                damping);

        this->constraints().push_back(std::move(constraint));
    }
}

void cXPBDTriangleMesh::computeNormals(void)
{
    N_.resize(num_faces,3);
    igl::per_face_normals(pdes_,F_,N_);
}
