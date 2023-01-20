#include "collision.h"

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
    }
}


// brute force broad phase
__global__ void simple_broad_phase(aabb* lhs, aabb* rhs, int* col, int ntris)
{
    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < ntris)
    {

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
}