#ifndef CUDAISFUN_SIMULATOR_H
#define CUDAISFUN_SIMULATOR_H

#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;


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
    void createTetrahedralMesh(float scale);

    // this function computes the normals in a parallel fashion
    void computeNormals(void);

    // host data pointers
    float* h_x; float* h_xlast; float* h_x0; float* h_xdot; float* h_m; int* h_ti; int* h_ei;

    // device data pointers
    int d_nx; float* d_x; float* d_xlast; float* d_x0; float* d_Nlast; float* d_N;
    float* d_xdot; int* d_ti; float* d_m; int d_ne; int* d_ei;

    // vis stuff
    int* h_fi; int* d_fi;

    // numbers
    int nvertices;
    int nelements;
    int nfaces;
    int nedges;

    // elements in eigen form
    Matrix<float,Dynamic,Dynamic,RowMajor> x;
    Matrix<float,Dynamic,Dynamic,RowMajor> xlast;
    Matrix<float,Dynamic,Dynamic,RowMajor> xdot;
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
    float h_young_modulus = 1000000; float h_poisson_ratio = 0.49; float h_mu; float h_lambda;
    float h_alpha; float h_beta; float* h_DmInv; float* h_v0;

    // define device variables
    float d_young_modulus; float d_poisson_ratio; float d_mu; float d_lambda;
    float d_alpha; float d_beta; float* d_DmInv; float* d_v0; float* d_lagrange;

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
    float h_alpha; float h_beta; float* h_v0;

    // define device variablesfloat d_young_modulus; float d_poisson_ratio; float d_mu; float d_lambda;
    float d_alpha; float d_beta; float* d_v0; float* d_lagrange;
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
    float h_alpha; float h_beta; float* h_e0;

    // define device variablesfloat d_young_modulus; float d_poisson_ratio; float d_mu; float d_lambda;
    float d_alpha; float d_beta; float* d_e0; float* d_lagrange;
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
    int* d_i; float* d_lagrange;
};


class simulator
{

public:

    simulator()
    {
        object = new meshObject();

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
    Matrix<float,Dynamic,Dynamic,RowMajor> positions();
    MatrixXi triangles(){return object->F;};

    // this funciton allocates memory to the fpu
    void allocMemory(void);

    // this function initializes a neoHookean constraint
    void initializeNeoHookeanConstraint(float alpha, float beta);

    // this function initializes a volume constraint
    void initializeVolumeConstraint(float alpha, float beta);

    // this funciton innitializes a edge length constraint
    void initializeEdgeConstraint(float alpha, float beta);

    // this function initializes a fixed point constraint
    void initializeFixedPointConstraint(void);

    // this function computes a neohookean constraint
    void computeNeohookeanConstraint(float* d_p);

    // this function computes the volume constraint
    void computeVolumeConstraint(float* d_p);

    // this function computes an edge constraint
    void computeEdgeConstraint(float* d_p);

    // this function computes a fixed point constraint
    void computeFixedPointConstraint(float* d_p);

    // this function computes collision constraint solving
    // including collision detection
    bool computeCollisionConstraints(float* pos, float* pos_last, float* proxy);

    // this function updates the dynamics of the mesh
    void updateDynamics(float* pos, float* pos_last, float* proxy, float* f, int i ,double dt);

    // this funciton frees the GPU
    void freeGPU(void);

public:

    // the mesh object
    meshObject* object;

    // neohookean constraint
    NeohookeanConstraint* nhc;

    // fixed point constraint
    FixedPointConstraint* fpc;

    // edge constraint
    EdgeConstraint* ec;

    // volume constraint
    VolumeConstraint* vc;

    // timestep (set 0.001 as default)
    float h_dt = 0.001;
    float d_dt = 0.001;

    // sim variable
    float* d_p;


};

#endif //CUDAISFUN_SIMULATOR_H
