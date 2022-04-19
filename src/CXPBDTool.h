//
// Created by aldof on 4/11/2022.
//

#include "chai3d.h"
#include "CXPBDConstraint.h"
#include "CXPBDAABB.h"
#include <Eigen/Core>

using namespace chai3d;

#ifndef CXPBD_CXPBDTOOL_H
#define CXPBD_CXPBDTOOL_H

class cXPBDTool : public cMesh
{

public:
    using positions_type   = Eigen::MatrixXd;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::MatrixX3d;
    using faces_type       = Eigen::MatrixXi;
    using elements_type    = Eigen::MatrixXi;
    using constraints_type = std::vector<std::unique_ptr<cXPBDConstraint>>;
    using scalar_type      = typename cXPBDConstraint::scalar_type;

    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of cMesh.
    cXPBDTool(cWorld* a_parentWorld) : cMesh()
    {
        m_tool = new cGenericTool(a_parentWorld);
    };

    //! Destructor of cMesh.
    ~cXPBDTool() {};

    //--------------------------------------------------------------------------
    // PUBLIC METHODS - GENERAL
    //--------------------------------------------------------------------------

    positions_type& positions() { return p_; }
    positions_type& positions_last() {return plast_;};
    int& numVerts(){return num_vertices;}
    faces_type& faces() { return F_; }
    int& numFaces(){return num_faces;}
    velocities_type& velocity() { return v_; }
    masses_type& mass() { return m_; }
    std::shared_ptr<AABBNode>& bb(){return bb_;}

    // This method connects the mesh to chai3d
    void connectToChai3d(void);

    // This method updates the chai3d mesh
    void updateChai3d(void);

    // This method computes the centroid
    Eigen::Vector3d computeCentroid();

    // This method translates the object
    void setLocalPos(Eigen::Vector3d a_pos);

    // This method sets the vertices for the object
    void setVertices(positions_type a_vertices)
    {
        p0_ = a_vertices;
        p_ = a_vertices;
        num_vertices = a_vertices.rows();
    };

    // This method sets the last set of vertices
    void setVerticesLast(void)
    {
        plast_ = p_;
    }

    // This method scales the object
    void scaleObject(double a_scale);

    // This method sets the velocities of the object
    void setVelocities(velocities_type a_velocity)
    {
        v_ = a_velocity;
    };

    // This method sets the outside faces for this object
    void setFaces(faces_type a_faces)
    {
        F_ = a_faces;
        num_faces = a_faces.rows();
    };

    // This method sets the mass of the deformable object
    void setMass(Eigen::VectorXd a_mass)
    {
        m_ = a_mass;
    }

    // This method builds a AABB boundary box
    void buildAABBBoundaryBox(void)
    {
        bb_ = buildAABB(plast_,p_,F_);
    }

    std::tuple<Eigen::MatrixX3d, // old position
            Eigen::MatrixX3d> // new position
    getPBDInfo()
    {
        int nverts = num_vertices;

        Eigen::MatrixXd oldPos(nverts, 3);
        Eigen::MatrixXd newPos(nverts, 3);

        oldPos = positions();
        newPos = positions();

        return std::make_tuple(oldPos, newPos);
    }

    void updatePos()
    {
        Eigen::Vector3d pos = m_tool->getDeviceLocalPos().eigen();
        setLocalPos(pos);
    }

public:
    cGenericTool* m_tool;

protected:
    positions_type const& p0() const { return p0_; }

private:
    positions_type p0_;            ///< Rest positions
    positions_type p_;             ///< Positions
    positions_type plast_;         ///< Last positions
    int num_vertices;              ///< Number of vertices
    faces_type F_;                 ///< Faces
    int num_faces;                 ///< Number of faces
    masses_type m_;                ///< Per-vertex mass_when_unfixed
    velocities_type v_;            ///< Per-vertex velocity
    std::shared_ptr<AABBNode> bb_; ///< Boundary box for this object



};

#endif //CXPBD_CXPBDTOOL_H
