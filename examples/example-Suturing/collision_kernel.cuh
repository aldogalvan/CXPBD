#include "aabb.cuh"

// brute force bounding box
__global__ void simple_bounding_box(aabb* aabbs, int* f, float* start, float* end, int ntris)
{

    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < ntris) {

        for (int i = 0; i < 3 ; i++)
        {
            aabbs[globalIdx].lower[i] = 999999;
            aabbs[globalIdx].upper[i] = -999999;
        }

        for (int k = 0 ; k < 3; k++)
        {
            int pt = f[3*globalIdx + k];
            for (int j = 0 ; j < 3; j++)
            {
                aabbs[globalIdx].lower[j] = min(aabbs[globalIdx].lower[j], start[3*pt + j]);
                aabbs[globalIdx].upper[j] = max(aabbs[globalIdx].upper[j], start[3*pt + j]);
                aabbs[globalIdx].lower[j] = min(aabbs[globalIdx].lower[j], end[3*pt + j]);
                aabbs[globalIdx].upper[j] = max(aabbs[globalIdx].upper[j], end[3*pt + j]);

            }
        }
        printf("%i, %f , %f , %f , %f , %f , %f \n",globalIdx , aabbs[globalIdx].upper[0],aabbs[globalIdx].upper[1],aabbs[globalIdx].upper[2],aabbs[globalIdx].lower[0],aabbs[globalIdx].lower[1],aabbs[globalIdx].lower[2] );
    }
}


// brute force broad phase
__global__ void simple_broad_phase(aabb* lhs, aabb* rhs, int* col, int ntris)
{
    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    //printf("RHS : %f , %f , %f , %f , %f , %f \n",rhs->lower.x , rhs->lower.y  , rhs->lower.z , rhs->upper.x , rhs->upper.y  , rhs->upper.z  );
    //printf("LHS : %f , %f , %f , %f , %f , %f \n",lhs[globalIdx].lower.x , lhs[globalIdx].lower.y  , lhs[globalIdx].lower.z , lhs[globalIdx].upper.x , lhs[globalIdx].upper.y  , lhs[globalIdx].upper.z  );

    if (globalIdx < ntris)
    {

        // printf("%f , %f , %f , %f , %f , %f \n",lhs[globalIdx].lower.x , lhs[globalIdx].lower.y  , lhs[globalIdx].lower.z , lhs[globalIdx].upper.x , lhs[globalIdx].upper.y  , lhs[globalIdx].upper.z  );
        for (int i = 0 ; i < 3; i++)
        {
            if (lhs[globalIdx].upper[i] < rhs->lower[i] || rhs->upper[i] < lhs[globalIdx].lower[i]) {
                col[globalIdx] = false;
                return;
            }
        }

        // if false then last is true
        col[globalIdx] = true;
    }

    //printf("%i - collision\n",globalIdx);
    //printf("RHS : %f , %f , %f , %f , %f , %f \n",rhs->lower.x , rhs->lower.y  , rhs->lower.z , rhs->upper.x , rhs->upper.y  , rhs->upper.z  );
    //printf("LHS : %f , %f , %f , %f , %f , %f \n",lhs[globalIdx].lower.x , lhs[globalIdx].lower.y  , lhs[globalIdx].lower.z , lhs[globalIdx].upper.x , lhs[globalIdx].upper.y  , lhs[globalIdx].upper.z  );
}