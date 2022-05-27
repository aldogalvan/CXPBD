//
// Created by aldof on 4/2/2022.
//

#include "chai3d.h"
#include "CXPBDConstraint.h"
#include "CXPBDAABB.h"
#include <set>

#ifndef CXPBD_CXPBDDEFORMABLEOBJECT_H
#define CXPBD_CXPBDDEFORMABLEOBJECT_H

using namespace std;
using namespace chai3d;

class cXPBDDeformableMesh : public cMesh
{

public:
    using positions_type   = Eigen::MatrixXd;
    using index_type       = int;
    using masses_type      = Eigen::VectorXd;
    using velocities_type  = Eigen::MatrixXd;
    using faces_type       = Eigen::MatrixXi;
    using elements_type    = Eigen::MatrixXi;
    using constraints_type = std::vector<std::unique_ptr<cXPBDConstraint>>;
    using scalar_type      = typename cXPBDConstraint::scalar_type;

    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    //! Constructor of cMesh.
    cXPBDDeformableMesh(){};

    //! Destructor of cMesh.
    ~cXPBDDeformableMesh(){};

    //--------------------------------------------------------------------------
    // PUBLIC METHODS - GENERAL
    //--------------------------------------------------------------------------

    positions_type& positions() { return p_; }
    positions_type& positions_last() {return plast_;};
    index_type& numVerts(){return num_vertices;}
    faces_type& faces() { return F_; }
    index_type& numFaces(){return num_faces;}
    elements_type& edges() { return E_; }
    index_type& numEdges(){return num_edges;}
    elements_type& tetrahedra() { return T_; }
    index_type& numTetrahedra(){return num_tetrahedra;}
    constraints_type& constraints() { return constraints_; }
    velocities_type& velocity() { return v_; }
    masses_type& mass() { return m_; }
    std::vector<bool>& fixed() { return fixed_; }
    std::shared_ptr<AABBNode>& bb(){return bb_;}

    // This method connects the mesh to chai3d
    void connectToChai3d(void);

    // This method updates the chai3d mesh
    void updateChai3d(Eigen::MatrixXd& a_pos);

    // Set frozen vertices
    void isFixed(vector<bool> a_fixed)
    {
        fixed_ = a_fixed;
    }

    // This method sets the ounding box
    void buildAABBBoundaryBox(positions_type pdes_)
    {
        bb_ = buildAABB(p_,pdes_,F_);
    }

    // This method computes the centroid
    Eigen::Vector3d computeCentroid();

    // This method translates the object
    void setLocalPos(Eigen::Vector3d a_pos);

    // This method scales the object
    void scaleObject(double a_scale);

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

    // This method sets the tetrahedra for the object
    void setTetrahedra(elements_type a_tetrahedra)
    {
        T_ = a_tetrahedra;
        num_tetrahedra = a_tetrahedra.rows();
    };

    // This method sets the edges for this object
    void setEdges(elements_type a_edges)
    {
        E_ = a_edges;
        num_edges = a_edges.rows();
    };

    // This method sets the mass of the deformable object
    void setMass(Eigen::VectorXd a_mass)
    {
        m_ = a_mass;
    }

    // This method sets the inside and outside vertices
    void setInsideOutside(set<index_type> in, set<index_type> out)
    {
        inside = in;
        in_size = in.size();
        outside = out;
        out_size = out.size();

    }

    // This method connects the tetrahedral mesh to the outside mesh
    void connectToSurface(void);

    // This method projects the tetrahedral mesh to the outside mesh
    void projectToSurface(void);

    // This method constrains the edge lengths
    void constrain_edge_lengths(scalar_type const compliance = 0.0);

    // This method constraints the volume
    void constrain_tetrahedron_volumes(scalar_type const compliance = 0.0);

    // This method constrains the hinges
    void constrain_hinge_bending(scalar_type const compliance = 0.0);

    // This method constrains the elastic potential
    void constrain_green_strain_elastic_potential(
            scalar_type young_modulus,
            scalar_type poisson_ratio,
            scalar_type const compliance = 0.0);

    // This method constrains the elasticity potential
    void constrain_neohookean_elasticity_potential(
            scalar_type young_modulus,
            scalar_type poisson_ratio,
            scalar_type const compliance = 0.0);


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

protected:

    positions_type const& p0() const { return p0_; }

private:

    positions_type p0_;            ///< Rest positions
    positions_type plast_;         ///< Last set of positions
    positions_type p_;             ///< Positions
    index_type num_vertices;       ///< Number of vertices
    faces_type F_;                 ///< Faces
    index_type num_faces;          ///< Number of faces
    elements_type E_;              ///< Elements
    index_type num_edges;          ///< Number of edges
    elements_type T_;              ///< Tetrahedra
    index_type num_tetrahedra;     ///< Number of tetrahedra
    elements_type H_;              ///< Hinges
    index_type num_hinges;         ///< Number of hinges
    masses_type m_;                ///< Per-vertex mass_when_unfixed
    velocities_type v_;            ///< Per-vertex velocity
    constraints_type constraints_; ///< PBD constraints
    std::vector<bool> fixed_;      ///< Flags fixed positions
    std::shared_ptr<AABBNode> bb_; ///< Boundary box for this object
    set<index_type> inside;        ///< Set of inside vertices
    index_type in_size;            ///< Size of inside set
    set<index_type> outside;       ///< Set of outside vertice
    index_type out_size;           ///< Size of outside set


};

#endif //CXPBD_CXPBDDEFORMABLEOBJECT_H
