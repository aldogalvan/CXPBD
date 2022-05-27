//
// Created by aldof on 4/2/2022.
//

#include "CXPBDDeformableObject.h"
#include "CXPBDEdgeLengthConstraint.h"
#include "CXPBDVolumeConstraint.h"
#include "CXPBDNeoHookeanConstraint.h"
#include "CXPBDBendingConstraint.h"

void cXPBDDeformableMesh::connectToChai3d(void)
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

void cXPBDDeformableMesh::updateChai3d(Eigen::MatrixXd& a_pos)
{
    for (int i = 0; i < num_vertices; i++)
    {
        Eigen::Vector3d vert = a_pos.row(i);
        m_vertices->setLocalPos(i,vert.x(),vert.y(),vert.z());
    }
}

Eigen::Vector3d cXPBDDeformableMesh::computeCentroid(void)
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

void cXPBDDeformableMesh::setLocalPos(Eigen::Vector3d a_pos)
{
    Eigen::Vector3d centroid = computeCentroid();
    Eigen::Vector3d trans = a_pos - centroid;

    for (int i = 0 ; i < num_vertices ; i++)
    {
        p0_.row(i) += trans;
        p_.row(i) += trans;
    }

}

void cXPBDDeformableMesh::scaleObject(double a_scale)
{
    p0_ *= a_scale;
    plast_ *= a_scale;
    p_ *= a_scale;
}

void cXPBDDeformableMesh::constrain_edge_lengths(scalar_type const compliance)
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
                compliance);

        this->constraints().push_back(std::move(constraint));
    }
}

void cXPBDDeformableMesh::constrain_tetrahedron_volumes(scalar_type const compliance)
{
    auto const& positions = this->p0();
    auto const& elements  = this->tetrahedra();

    for (auto i = 0u; i < elements.rows(); ++i)
    {
        auto const element = elements.row(i);
        auto constraint    = std::make_unique<cXPBDVolumeConstraint>(
                std::initializer_list<std::uint32_t>{
                        static_cast<std::uint32_t>(element(0)),
                        static_cast<std::uint32_t>(element(1)),
                        static_cast<std::uint32_t>(element(2)),
                        static_cast<std::uint32_t>(element(3))},
                positions,
                compliance);

        this->constraints().push_back(std::move(constraint));
    }
}

void cXPBDDeformableMesh::constrain_hinge_bending(const scalar_type compliance)
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
                compliance);

        this->constraints().push_back(std::move(constraint));
    }
}

void cXPBDDeformableMesh::constrain_neohookean_elasticity_potential(
        scalar_type young_modulus,
        scalar_type poisson_ratio,
        scalar_type const compliance)
{
    auto const& positions = this->p0();
    auto const& elements  = this->tetrahedra();

    for (auto i = 0u; i < elements.rows(); ++i)
    {
        auto const element = elements.row(i);
        auto constraint    = std::make_unique<cXPBDNeoHookeanConstraint>(
                std::initializer_list<std::uint32_t>{
                        static_cast<std::uint32_t>(element(0)),
                        static_cast<std::uint32_t>(element(1)),
                        static_cast<std::uint32_t>(element(2)),
                        static_cast<std::uint32_t>(element(3))},
                positions,
                young_modulus,
                poisson_ratio,
                compliance);

        this->constraints().push_back(std::move(constraint));
    }
}