//
// Created by aldof on 4/11/2022.
//

#include "CXPBDTool.h"

void cXPBDTool::connectToChai3d(void)
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

void cXPBDTool::updateChai3d()
{
    for (int i = 0; i < num_vertices; i++)
    {
        Eigen::Vector3d vert = p_.row(i);
        m_vertices->setLocalPos(i,vert.x(),vert.y(),vert.z());
    }
}

Eigen::Vector3d cXPBDTool::computeCentroid(void)
{
    Eigen::Vector3d sum;
    sum.setZero();

    for (int i = 0 ; i < num_vertices ; i++)
    {
        sum += p0_.row(i);
    }

    sum /= num_vertices;
}

void cXPBDTool::scaleObject(double a_scale)
{
    p0_ *= a_scale;
    p_ *= a_scale;
}


void cXPBDTool::setLocalPos(Eigen::Vector3d a_pos)
{
    Eigen::Vector3d centroid = computeCentroid();
    Eigen::Vector3d trans = a_pos - centroid;

    for (int i = 0 ; i < num_vertices ; i++)
    {
        p_.row(i) += trans;
    }

}
