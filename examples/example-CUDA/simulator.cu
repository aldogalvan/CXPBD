#include "simulator_kernel.cuh"
#include "simulator.h"
#include "collision.h"


// returns the position data of the faces
Matrix<double,Dynamic,Dynamic,RowMajor> simulator::positions()
{

    double* temp;
    temp = (double*)malloc(object->nvertices*3*sizeof (double));
    cudaMemcpy(temp,object->d_x,3*object->nvertices*sizeof(double),cudaMemcpyDeviceToHost);
    Matrix<double,Dynamic,Dynamic,RowMajor> ret(object->nvertices,3);
    for (int i = 0; i < object->nvertices; i++)
    {
        ret(i,0) = temp[3*i+0];
        ret(i,1) = temp[3*i+1];
        ret(i,2) = temp[3*i+2];
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
    cudaMalloc((void**)&this->d_x,this->nvertices*3*sizeof(double));
    cudaMalloc((void**)&this->d_xlast,this->nvertices*3*sizeof(double));
    cudaMalloc((void**)&this->d_x0,this->nvertices*3*sizeof(double));
    cudaMalloc((void**)&this->d_xdot,this->nvertices*3*sizeof(double));
    cudaMalloc((void**)&this->d_N,this->nfaces*3*sizeof(double));
    cudaMalloc((void**)&this->d_Nlast,this->nfaces*3*sizeof(double));
    cudaMalloc((void**)&this->d_m,this->nvertices*sizeof(double));


    cudaMemcpy(this->d_ti,this->h_ti,4*this->nelements*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_ei,this->h_ei,2*this->nedges*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_fi,this->h_fi,3*this->nfaces*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_x,this->h_x,3*this->nvertices*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_xlast,this->h_x,3*this->nvertices*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_x0,this->h_x0,3*this->nvertices*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_xdot,this->h_xdot,3*this->nvertices*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_m,this->h_m,this->nvertices*sizeof(double),cudaMemcpyHostToDevice);
}

void NeohookeanConstraint::transferToGPU()
{
    // defines constants
    d_mu = h_mu; d_lambda = h_lambda; d_alpha = h_alpha; d_beta = h_beta;

    //allocate memory
    cudaMalloc((void**)&this->d_DmInv,9*this->h_nconstraints*sizeof(double));
    cudaMalloc((void**)&this->d_v0,this->h_nconstraints*sizeof(double));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(double));


    // copy to device
    cudaMemcpy(this->d_DmInv,this->h_DmInv,9*this->h_nconstraints*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_v0,this->h_v0,this->h_nconstraints*sizeof(double),cudaMemcpyHostToDevice);
}

void VolumeConstraint::transferToGPU()
{
    // defines constants
    d_alpha = h_alpha; d_beta = h_beta;

    //allocate memory
    cudaMalloc((void**)&this->d_v0,this->h_nconstraints*sizeof(double));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(double));
    cudaMalloc((void**)&this->d_graph,this->maxcolor*this->maxelem*sizeof(int));

    // copy to device
    cudaMemcpy(this->d_v0,this->h_v0,this->h_nconstraints*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_graph,this->h_graph,this->maxcolor*this->maxelem*sizeof(int),cudaMemcpyHostToDevice);
}

void EdgeConstraint::transferToGPU()
{

    // defines constants
    d_alpha = h_alpha; d_beta = h_beta;

    //allocate memory
    cudaMalloc((void**)&this->d_e0,this->h_nconstraints*sizeof(double));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(double));
    cudaMalloc((void**)&this->d_graph,this->maxcolor*this->maxelem*sizeof(int));

    // copy to device
    cudaMemcpy(this->d_e0,this->h_e0,this->h_nconstraints*sizeof(double),cudaMemcpyHostToDevice);
    cudaMemcpy(this->d_graph,this->h_graph,this->maxcolor*this->maxelem*sizeof(int),cudaMemcpyHostToDevice);
}

void FixedPointConstraint::transferToGPU()
{
    //allocate memory
    cudaMalloc((void**)&this->d_i,this->h_nconstraints*sizeof(int));
    cudaMalloc((void**)&this->d_lagrange,this->h_nconstraints*sizeof(double));

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
    cudaMalloc((void**)&d_p,object->nvertices*3*sizeof(double));;
}

void simulator::computeNeohookeanConstraint(double* d_p)
{

    neohookeanConstraint<<<1,1024>>>(d_dt,nhc->h_nconstraints,nhc->d_mu,nhc->d_lambda,
                                  nhc->d_alpha,nhc->d_beta,object->d_ti, d_p,object->d_x0,
                                  nhc->d_DmInv,nhc->d_v0, object->d_m, nhc->d_lagrange);
    cudaDeviceSynchronize();

}

void simulator::computeEdgeConstraint(double* d_p)
{
    for (int c = 0; c < ec->maxcolor; c++)
    {
        edgeConstraint<<<4, 1024>>>(d_dt, ec->h_nconstraints, ec->d_alpha, ec->d_beta, object->d_ei, d_p, object->d_x0,
                                    ec->d_e0, object->d_m, ec->d_lagrange, ec->d_graph, c, ec->maxelem);
    }
    cudaDeviceSynchronize();
}

void simulator::computeVolumeConstraint(double* d_p)
{

    for (int c = 0; c < vc->maxcolor; c++)
    {
        volumeConstraint<<<4, 1024>>>(d_dt, vc->h_nconstraints, vc->d_alpha, vc->d_beta, object->d_ti, d_p,
                                      object->d_x0, vc->d_v0, object->d_m, vc->d_lagrange, vc->d_graph, c, vc->maxelem);

    }
    cudaDeviceSynchronize();
}

void simulator::computeFixedPointConstraint(double* d_p)
{

    fixedPointConstraint<<<4,1024>>>(fpc->h_nconstraints,fpc->d_i,d_p,object->d_x0);
    cudaDeviceSynchronize();
}

bool simulator::computeCollisions(vector<ColInfo*>& col_info)
{
    // first find any collisions
    return findCollisions(tool,object,col_info);
}

void simulator::computeCollisionConstraints(vector<ColInfo*>& col_info)
{
    // penalty method
    auto collision = col_info[0];
    auto normal = collision->normal[1];
    Vector3f pos = Vector3f(tool->pos[0],tool->pos[1],tool->pos[2]);
    double k = tool->k_proxy;
    auto triangle = collision->end;
    auto a = triangle[0];
    auto b = triangle[1];
    auto c = triangle[2];
    Vector3d ret = closestPointOnTriangle(triangle,pos.cast<double>());

    // cout << (pos.cast<double>() - ret).norm() << endl;
    double distance = (pos.cast<double>() - ret).norm();
    if ( distance < tool->tool_radius)
    {
        if ((pos.cast<double>() - a).dot(normal) > 0)
        {
            auto ret_force = distance*normal*tool->k_proxy;
            tool->force[0] = ret_force[0];
            tool->force[1] = ret_force[1];
            tool->force[2] = ret_force[2];

            tool->force_mesh[0] = -ret_force[0];
            tool->force_mesh[1] = -ret_force[1];
            tool->force_mesh[2] = -ret_force[2];

            tool->index = collision->collision_index_;

        }
        else
        {
            tool->force[0] = 0;
            tool->force[1] = 0;
            tool->force[2] = 0;
        }
    }
    else
    {
        tool->force[0] = 0;
        tool->force[1] = 0;
        tool->force[2] = 0;
    }


    /*
    // error measure
    double ep = 1e-6;

    // do the haptics!
    // you got this Aldo!

    // implicit collison handling
    //*************************************************

    // do case where the proxy collides with the object

    // check if it is an inside or outside collision


    Vector3d normal_col = normal_start + t_*(normal_end - normal_start);
    Vector3d q_col = q_start + t_*(q_end - q_start);
    Vector3d a_col = a_start + t_*(a_end - a_start);
    Vector3d b_col = b_start + t_*(b_end - b_start);
    Vector3d c_col = c_start + t_*(c_end - c_start);

    // barycentric coordinates for time of collision
    double barycentric_col[4];
    Barycentric(barycentric_col,a_col,b_col,c_col,q_col);

    // the desired position
    Vector3d proxy_temp = (barycentric_col[0]*a_end/barycentric_col[3] + barycentric_col[1]*b_end/barycentric_col[3]
                           + barycentric_col[2]*c_end/barycentric_col[3]);

    // assert that the distance between the proxy > toolRadius
    Vector3d mat[3];
    mat[0] = a_end;
    mat[1] = b_end;
    mat[2] = c_end;

    Vector3d pt = closestPointOnTriangle(mat,proxy_temp);
    std::cout << (proxy_temp - pt).norm() << std::endl;

    // move it a bit in the normal direction
    proxy_temp += ep*normal_end;
    proxy[0] = (double)proxy_temp[0]; proxy[1] = (double)proxy_temp[1]; proxy[2] = (double)proxy_temp[2];

    // do case where the object collides with the proxy
    //*************************************************

    // explicit collision handling
    //*************************************************
    // do case where the proxy collides with the object

    // do case where the object collides with the proxy
    //*************************************************
    ret = 1;
*/
 }

void simulator::updateDynamics(double dt)
{
    int idx = 100;
    Vector3i i;
    if (idx != -1)
        i = Vector3i(object->h_fi[3*idx+0],object->h_fi[3*idx+1],object->h_fi[3*idx+2]);
    Vector3d f(0,0,0);

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
    object->computeNormals();

    vector<ColInfo*> col_info;

    if(computeCollisions(col_info))
    {
        computeCollisionConstraints(col_info);
        tool->proxy_pos[0] = tool->pos[0]; tool->proxy_pos[1] = tool->pos[1]; tool->proxy_pos[2] = tool->pos[2];
    }
    else
    {
        // check if there was a collision
        // if no then the proxy updates to current position and force is zero
        tool->proxy_pos[0] = tool->pos[0]; tool->proxy_pos[1] = tool->pos[1]; tool->proxy_pos[2] = tool->pos[2];

        tool->force[0] = 0;
        tool->force[1] = 0;
        tool->force[2] = 0;
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