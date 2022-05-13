//
// Created by aldof on 5/5/2022.
//

#include "CXPBDToolMesh.h"


void cXPBDToolMesh::connectToChai3d(void)
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

void cXPBDToolMesh::updateChai3d()
{
    for (auto i = 0u; i < num_vertices; i++)
    {
        Eigen::Vector3d vert = p_.row(i);
        //Eigen::Vector3d vert = pnext_.row(i);
        m_vertices->setLocalPos(i,vert.x(),vert.y(),vert.z());
    }
}

void cXPBDToolMesh::computeCentroid(void)
{
    Eigen::Vector3d sum;
    sum.setZero();

    for (int i = 0 ; i < num_vertices ; i++)
    {
        sum += p0_.row(i);
    }

    sum /= num_vertices;

    pos = sum;
}

void cXPBDToolMesh::scaleObject(double a_scale)
{
    p0_ *= a_scale;
    p_ *= a_scale;
}


void cXPBDToolMesh::setLocalPos(Eigen::Vector3d a_pos)
{

    Eigen::Vector3d trans = a_pos - pos;

    for (int i = 0 ; i < num_vertices ; i++)
    {
        p_.row(i) += trans;
    }
    pos = a_pos;

}
