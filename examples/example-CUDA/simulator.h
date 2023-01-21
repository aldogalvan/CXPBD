#ifndef CUDAISFUN_SIMULATOR_H
#define CUDAISFUN_SIMULATOR_H

#include <vector>
#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

struct ColInfo;

struct toolObject
{
    toolObject(double tool_radius_, double k_proxy_, double b_proxy_) :
        tool_radius(tool_radius_), k_proxy(k_proxy_), b_proxy(b_proxy_)
    {
        force[0] = 0; force[1] = 0; force[2] = 0;
    }

    ~toolObject()
    {

    }

    double pos[3];
    double last_pos[3];
    double proxy_pos[3];
    double force[3];
    double tool_radius;
    double k_proxy;
    double b_proxy;

    // force and index
    double force_mesh[3];
    int index;
};

struct meshObject
{
    meshObject()
    {
        createTetrahedralMesh(1000.0);
        transferToGPU();
        computeNormals();

        cout << "Object Constructed" << endl;
    }

    ~meshObject()
    {
        freeGPU();
    }

    // this function transfers the data to the GPU
    // this function transfers data to the GPU
    void transferToGPU(void);

    // this function frees the GPU
    void freeGPU(void);

    // this function creates the mesh object
    void createTetrahedralMesh(double scale);

    // this function computes the normals in a parallel fashion
    void computeNormals(void);

    // host data pointers
    double* h_x; double* h_xlast; double* h_x0; double* h_xdot; double* h_m; int* h_ti; int* h_ei;

    // device data pointers
    int d_nx; double* d_x; double* d_xlast; double* d_x0; double* d_Nlast; double* d_N;
    double* d_xdot; int* d_ti; double* d_m; int d_ne; int* d_ei;

    // vis stuff
    int* h_fi; int* d_fi;

    // numbers
    int nvertices;
    int nelements;
    int nfaces;
    int nedges;

    // elements in eigen form
    Matrix<double,Dynamic,Dynamic,RowMajor> x;
    Matrix<double,Dynamic,Dynamic,RowMajor> xlast;
    Matrix<double,Dynamic,Dynamic,RowMajor> xdot;
    Matrix<int,Dynamic,Dynamic,RowMajor> F;
    Matrix<int,Dynamic,Dynamic,RowMajor> T;
    Matrix<int,Dynamic,Dynamic,RowMajor> E;

};

struct NeohookeanConstraint
{
    // this function transfers all contents to the GPU
    void transferToGPU(void);

    // this function computes a graph for this constraint
    void computeGraph(MatrixXi E);

    // number of constraints (elements)
    int h_nconstraints;
    int d_nconstraints;

    // graph
    int* h_graph;
    int* d_graph;

    // values for neohookean constraint
    double h_young_modulus = 1000000; double h_poisson_ratio = 0.49; double h_mu; double h_lambda;
    double h_alpha; double h_beta; double* h_DmInv; double* h_v0;

    // define device variables
    double d_young_modulus; double d_poisson_ratio; double d_mu; double d_lambda;
    double d_alpha; double d_beta; double* d_DmInv; double* d_v0; double* d_lagrange;

};

struct VolumeConstraint
{
    // this function transfers all contents to the GPU
    void transferToGPU(void);

    // this function computes a graph for this constraint
    void computeGraph(MatrixXi E);

    // number of constraints (elements)
    int h_nconstraints;
    int d_nconstraints;

    // graph
    int* h_graph;
    int* d_graph;
    int maxcolor;
    int maxelem;

    // values for neohookean constraint
    double h_alpha; double h_beta; double* h_v0;

    // define device variablesdouble d_young_modulus; double d_poisson_ratio; double d_mu; double d_lambda;
    double d_alpha; double d_beta; double* d_v0; double* d_lagrange;
};

struct EdgeConstraint
{
    // this function transfers all contents to the GPU
    void transferToGPU(void);

    // this function computes a graph for this constraint
    void computeGraph(MatrixXi E);

    // number of constraints (elements)
    int h_nconstraints;
    int d_nconstraints;

    // graph
    int* h_graph;
    int* d_graph;
    int maxcolor;
    int maxelem;

    // values for neohookean constraint
    double h_alpha; double h_beta; double* h_e0;

    // define device variablesdouble d_young_modulus; double d_poisson_ratio; double d_mu; double d_lambda;
    double d_alpha; double d_beta; double* d_e0; double* d_lagrange;
};

struct FixedPointConstraint
{
    // this function transfers all contents to the GPU
    void transferToGPU(void);

    // number of constraints (elements)
    int h_nconstraints;
    int d_nconstraints;

    // host variables
    int* h_i;

    // device variables
    int* d_i; double* d_lagrange;
};


class simulator
{

public:

    simulator(toolObject* a_tool)
    {
        object = new meshObject();
        tool = a_tool;

        initializeEdgeConstraint(0.99,0.0);
        initializeVolumeConstraint(0.99,0.0);
        initializeFixedPointConstraint();

        allocMemory();
        object->computeNormals();

    }

    ~simulator()
    {
        freeGPU();
        delete object;
    }

public:


    // this function returns a pointer containing current position data of
    //  the object
    Matrix<double,Dynamic,Dynamic,RowMajor> positions();
    MatrixXi triangles(){return object->F;};

    // this funciton allocates memory to the fpu
    void allocMemory(void);

    // this function initializes a neoHookean constraint
    void initializeNeoHookeanConstraint(double alpha, double beta);

    // this function initializes a volume constraint
    void initializeVolumeConstraint(double alpha, double beta);

    // this funciton innitializes a edge length constraint
    void initializeEdgeConstraint(double alpha, double beta);

    // this function initializes a fixed point constraint
    void initializeFixedPointConstraint(void);

    // this function computes a neohookean constraint
    void computeNeohookeanConstraint(double* d_p);

    // this function computes the volume constraint
    void computeVolumeConstraint(double* d_p);

    // this function computes an edge constraint
    void computeEdgeConstraint(double* d_p);

    // this function computes a fixed point constraint
    void computeFixedPointConstraint(double* d_p);

    // this function computes collision constraint solving
    // including collision detection
    bool computeCollisions(vector<ColInfo*>& col_info);

    // this function computes the collision constraints
    void computeCollisionConstraints(vector<ColInfo*>& col_info);

    // this function updates the dynamics of the mesh
    void updateDynamics(double dt);

    // this funciton frees the GPU
    void freeGPU(void);

public:

    // the mesh object
    meshObject* object;

    // the tool object
    toolObject* tool;

    // neohookean constraint
    NeohookeanConstraint* nhc;

    // fixed point constraint
    FixedPointConstraint* fpc;

    // edge constraint
    EdgeConstraint* ec;

    // volume constraint
    VolumeConstraint* vc;

    // timestep (set 0.001 as default)
    double h_dt = 0.001;
    double d_dt = 0.001;

    // sim variable
    double* d_p;

};

#endif //CUDAISFUN_SIMULATOR_H
