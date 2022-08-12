//
// Created by agalvan-admin on 8/2/22.
//

#include "CXPBDThread.h"
#include "constraints/CXPBDEdgeLengthConstraint.h"
#include "constraints/CXPBDDynamicFixedPointConstraint.h"

void cXPBDThread::setVerticesLast()
{
    plast_ = p_;
}

void cXPBDThread::connectToChai3d(void)
{

    for (int i = 0; i < num_vertices; i++)
    {
        Vector3d vert = p_.row(i);
        newVertex(vert.x(),vert.y(),vert.z());
    }

    for (int i = 0; i < num_edges; i++)
    {
        Vector2i edge = E_.row(i);
        m_edges.emplace_back(cEdge(this,edge(0),edge(1)));
    }

}

void cXPBDThread::updateChai3d(void)
{
    for (int i = 0; i < num_vertices; i++)
    {
        Eigen::Vector3d vert = p_.row(i);
        m_vertices->setLocalPos(i,vert.x(),vert.y(),vert.z());
    }
}

void cXPBDThread::setSuturePos(std::shared_ptr<Vector3d> a_pos)
{
    suture_pos = a_pos;
}

void cXPBDThread::constrain_dynamic_point(const cXPBDThread::scalar_type compliance,
                                          const cXPBDThread::scalar_type damping)
{

    auto constraint = std::make_unique<cXPBDDynamicFixedPointConstraint>(
            std::initializer_list<std::uint32_t>{
                    static_cast<std::uint32_t>(0)},
            suture_pos,
            compliance,
            damping);

    this->constraints().push_back(std::move(constraint));

}

void cXPBDThread::constrain_edge_bending(const cXPBDThread::scalar_type compliance,
                                         const cXPBDThread::scalar_type damping)
{

}

void cXPBDThread::constrain_edge_lengths(const cXPBDThread::scalar_type compliance,
                                         const cXPBDThread::scalar_type damping)
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