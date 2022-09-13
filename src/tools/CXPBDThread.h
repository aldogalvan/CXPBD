#include "chai3d.h"
#include <Eigen/Dense>
#include "../constraints/CXPBDConstraint.h"
#include "../collision/CXPBDAABB.h"
#include <set>

#ifndef CXPBD_CXPBDTHREAD_H
#define CXPBD_CXPBDTHREAD_H

using namespace chai3d;
using namespace Eigen;
using namespace std;

class cXPBDThread : public cMesh
{

public:

    using positions_type   = MatrixXd;
    using index_type       = int;
    using masses_type      = VectorXd;
    using velocities_type  = MatrixXd;
    using faces_type       = MatrixXi;
    using elements_type    = MatrixXi;
    using constraints_type = std::vector<std::unique_ptr<cXPBDConstraint>>;
    using scalar_type      = typename cXPBDConstraint::scalar_type;

    //--------------------------------------------------------------------------
    // CONSTRUCTOR & DESTRUCTOR:
    //--------------------------------------------------------------------------

public:

    cXPBDThread(int numSegments , double length , double radius )
    {
        // the number of edges in the thread
        num_edges = numSegments;
        // the number of points in the thread
        num_vertices = numSegments + 1;
        // the number of hinges in the thread
        num_hinges = numSegments - 1;
        // resize and initialize to zero
        p_.resize(num_vertices,3);
        p_.setZero(); pdes_.setZero(); plast_.setZero();
        // resize edge matrix
        E_.resize(num_edges,2);
        // resize hinge matrix
        H_.resize(num_hinges,3);
        double step = length / numSegments;

        for (int i = 0u; i < num_vertices; i++)
        {

            // draw thread along the y dimension
            p_.row(i) = i*step*Vector3d(0,1,0);
            // set the edge indices
            if (i < num_vertices - 1)
                E_.row(i) = Vector2i(i,i+1);
            // set the hinge indices
            if (i < num_vertices - 2)
                H_.row(i) = Vector3i(i,i+1,i+2);
        }

        pdes_ = p_; p0_ = p_; plast_ = p_;
        v_ = p_; v_.setZero();
        m_.resize(num_vertices);
        m_.setConstant(0.001);
    }

    ~cXPBDThread()
    {

    }

    //--------------------------------------------------------------------------
    // PUBLIC METHODS - GENERAL
    //--------------------------------------------------------------------------

public:


    positions_type& positions() { return p_; }
    positions_type& positions_last() {return plast_;}
    positions_type& positions_desired() {return pdes_;}
    positions_type& p0() { return p0_; }
    index_type& numVerts(){return num_vertices;}
    elements_type& edges() { return E_; }
    constraints_type& constraints() { return constraints_; }
    velocities_type& velocity() { return v_; }
    masses_type& mass() { return m_; }


public:


    //--------------------------------------------------------------------------
    // PUBLIC METHODS - GENERAL
    //--------------------------------------------------------------------------

    // This method sets the last set of vertices
    void setVerticesLast(void);

    // this method connects this object for chai3d visualiztion
    void connectToChai3d(void);

    // this method updates the chai3d object
    void updateChai3d(void);

    // this method sets a new suture position and ties it to the node of index 0
    void setSuturePos(std::shared_ptr<Vector3d> a_pos);

    // this method draws the thread from basic primitives from OpenGL/chai3d
    void drawThread(void);

    // this method creates a dynamic fixed point constraint
    void constrain_dynamic_point(scalar_type const compliance = 0.0 , scalar_type const damping = 0.0);

    // This method constrains the edge lengths
    void constrain_edge_lengths(scalar_type const compliance = 0.0 , scalar_type const damping = 0.0);

    // This method constrains the edge lengths
    void constrain_edge_bending(scalar_type const compliance = 0.0 , scalar_type const damping = 0.0);


private:

    int m_numSegments = 100;       ///< The number of segments in the thread
    double m_length = 1.00;        ///< Total length
    positions_type p0_;            ///< Rest positions
    positions_type plast_;         ///< Last set of positions
    positions_type p_;             ///< Positions
    positions_type pdes_;          ///< Desired positions
    index_type num_vertices;       ///< Number of vertices
    elements_type E_;              ///< Edges
    index_type num_edges;          ///< Number of edges
    elements_type H_;              ///< Hinges
    index_type num_hinges;         ///< Number of hinges
    masses_type m_;                ///< Per-vertex mass_when_unfixed
    velocities_type v_;            ///< Per-vertex velocity
    constraints_type constraints_; ///< PBD constraints
    scalar_type r_;                ///< Radius of the thread

private:

    std::shared_ptr<Vector3d> suture_pos;  ///< The position of the suture

};

#endif //CXPBD_CXPBDTHREAD_H
