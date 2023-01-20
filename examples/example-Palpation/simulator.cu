#include "simulator_kernel.cuh"
#include "simulator.h"
#include "collision.h"


// returns the position data of the faces
Matrix<float,Dynamic,Dynamic,RowMajor> simulator::positions()
{

    float* temp;
    temp = (float*)malloc(object->nvertices*3*sizeof (float));
    cudaMemcpy(temp,object->d_x,3*object->nvertices*sizeof(float),cudaMemcpyDeviceToHost);
    Matrix<float,Dynamic,Dynamic,RowMajor> ret(object->nvertices,3);

    for (int i = 0; i < object->nvertices; i++)
    {
        ret(i,0) = temp[3*i+0] / 1000;
        ret(i,1) = temp[3*i+1] / 1000;
        ret(i,2) = temp[3*i+2] / 1000;
    }

    return ret;
}

// ! Everything here will live in the GPU
void meshObject::transferToGPU(void)
{
    // define object variables
    this->d_ne = this->nelements;
    this->d_nx = this->nvertices;

    // allocate space for object variables
    cudaMalloc((void**)&this->d_ti,this->nelements*4*sizeof(int));
    cudaMalloc((void**)&this->d_ei,this->nedges*2*sizeof(int));
    cudaMalloc((void**)&this->d_fi,this->nfaces*3*sizeof(int));
    cudaMalloc((void**)&this->d_x,this->nvertices*3*sizeof(float));
    cudaMalloc((void**)&this->d_xlast,this->nvertices*3*sizeof(float));
    cudaMalloc((void**)&this->d_x0,this->nvertices*3*sizeof(float));
    cudaMalloc((void**)&this->d_xdot,this->nvertices*3*sizeof(float));
    cudaMalloc((void**)&this->d_N,this->nfaces*3*sizeof(float));
    cudaMalloc((void**)&this->d_Nlast,this->nfaces*3*sizeof(float));
    cudaMalloc((void**)&this->d_m,this->nvertices*sizeof(float));


    cudaMemcpy(this->d_ti,this->h_ti,4*this->nelements*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_ei,this->h_ei,2*this->nedges*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_fi,this->h_fi,3*this->nfaces*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_x,this->h_x,3*this->nvertices*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_xlast,this->h_x,3*this->nvertices*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_x0,this->h_x0,3*this->nvertices*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_xdot,this->h_xdot,3*this->nvertices*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_m,this->h_m,this->nvertices*sizeof(float),cudaMemcpyHostToDevice);
}

void NeohookeanConstraint::transferToGPU()
{
    // defines constants
    d_mu = h_mu; d_lambda = h_lambda; d_alpha = h_alpha; d_beta = h_beta;

    //allocate memory
    cudaMalloc((void**)&this->d_DmInv,9*this->h_nconstraints*sizeof(float));
    cudaMalloc((void**)&this->d_v0,this->h_nconstraints*sizeof(float));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(float));


    // copy to device
    cudaMemcpy(this->d_DmInv,this->h_DmInv,9*this->h_nconstraints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_v0,this->h_v0,this->h_nconstraints*sizeof(float),cudaMemcpyHostToDevice);
}

void VolumeConstraint::transferToGPU()
{
    // defines constants
    d_alpha = h_alpha; d_beta = h_beta;

    //allocate memory
    cudaMalloc((void**)&this->d_v0,this->h_nconstraints*sizeof(float));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(float));
    cudaMalloc((void**)&this->d_graph,this->maxcolor*this->maxelem*sizeof(int));

    // copy to device
    cudaMemcpy(this->d_v0,this->h_v0,this->h_nconstraints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_graph,this->h_graph,this->maxcolor*this->maxelem*sizeof(int),cudaMemcpyHostToDevice);
}

void EdgeConstraint::transferToGPU()
{

    // defines constants
    d_alpha = h_alpha; d_beta = h_beta;

    //allocate memory
    cudaMalloc((void**)&this->d_e0,this->h_nconstraints*sizeof(float));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(float));
    cudaMalloc((void**)&this->d_graph,this->maxcolor*this->maxelem*sizeof(int));

    // copy to device
    cudaMemcpy(this->d_e0,this->h_e0,this->h_nconstraints*sizeof(float),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_graph,this->h_graph,this->maxcolor*this->maxelem*sizeof(int),cudaMemcpyHostToDevice);
}

void FixedPointConstraint::transferToGPU()
{
    //allocate memory
    cudaMalloc((void**)&this->d_i,this->h_nconstraints*sizeof(int));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(float));

    // copy to device
    cudaMemcpy(this->d_i,this->h_i,this->h_nconstraints*sizeof(int),cudaMemcpyHostToDevice);
}


void meshObject::computeNormals(void)
{
    ::computeNormals<<<2,1024>>>(this->d_x,this->d_Nlast, this->d_N,this->d_fi,this->nfaces);
    cudaDeviceSynchronize();
}

void simulator::allocMemory(void)
{
    cudaMalloc((void**)&d_p,object->nvertices*3*sizeof(float));;
}

void simulator::computeNeohookeanConstraint(float* d_p)
{

    neohookeanConstraint<<<1,1024>>>(d_dt,nhc->h_nconstraints,nhc->d_mu,nhc->d_lambda,
                                  nhc->d_alpha,nhc->d_beta,object->d_ti, d_p,object->d_x0,
                                  nhc->d_DmInv,nhc->d_v0, object->d_m, nhc->d_lagrange);
    cudaDeviceSynchronize();

}

void simulator::computeEdgeConstraint(float* d_p)
{

    for (int c = 0; c < ec->maxcolor; c++)
    {
        edgeConstraint<<<4, 1024>>>(d_dt, ec->h_nconstraints, ec->d_alpha, ec->d_beta, object->d_ei, d_p, object->d_x0,
                                    ec->d_e0, object->d_m, ec->d_lagrange, ec->d_graph, c, ec->maxelem);
    }
    cudaDeviceSynchronize();
}

void simulator::computeVolumeConstraint(float* d_p)
{

    for (int c = 0; c < vc->maxcolor; c++)
    {
        volumeConstraint<<<4, 1024>>>(d_dt, vc->h_nconstraints, vc->d_alpha, vc->d_beta, object->d_ti, d_p,
                                      object->d_x0, vc->d_v0, object->d_m, vc->d_lagrange, vc->d_graph, c, vc->maxelem);

    }
    cudaDeviceSynchronize();
}

void simulator::computeFixedPointConstraint(float* d_p)
{

    fixedPointConstraint<<<4,1024>>>(fpc->h_nconstraints,fpc->d_i,d_p,object->d_x0);
    cudaDeviceSynchronize();
}

bool simulator::computeCollisionConstraints(float* pos, float* pos_last, float* proxy)
{
    // first find any collisions
    return findCollisions(pos,pos_last,proxy,object);
}

void simulator::updateDynamics(float* pos, float* pos_last, float* proxy, float* f, int i, double dt)
{

    // first we timestep using explicit euler
    explicitEuler<<<4,1024>>>(i, f, object->d_nx, d_dt,
            d_p, object->d_x, object->d_xdot, object->d_m );

    cudaDeviceSynchronize();

    // compute the constraints for the object
    computeVolumeConstraint(d_p);
    computeEdgeConstraint(d_p);
    computeFixedPointConstraint(d_p);

    // updates the solution
    updateSolution<<<4,1024>>>(d_dt, object->d_nx, object->d_xdot, object->d_x,
                               object->d_xlast,d_p);
    cudaDeviceSynchronize();

    // find collisions
    // TODO: Implement actual tree AABB, this is brute forcing but fine for now
    object->computeNormals();

    if(computeCollisionConstraints(pos,pos_last,proxy))
    {

    }
}


void meshObject::freeGPU(void)
{
    cudaFree(this->d_ti);
    cudaFree(this->d_ei);
    cudaFree(this->d_fi);
    cudaFree(this->d_x);
    cudaFree(this->d_xlast);
    cudaFree(this->d_x0);
    cudaFree(this->d_xdot);
    cudaFree(this->d_N);
    cudaFree(this->d_Nlast);
    cudaFree(this->d_m);

}

void simulator::freeGPU(void)
{
    cudaFree(d_p);
}