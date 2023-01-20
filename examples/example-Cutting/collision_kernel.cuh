#include "aabb.cuh"

// brute force bounding box
__global__ void simple_bounding_box(aabb** aabbs, int* f, float* start, float* end, int ntris)
{

    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < ntris) {

        aabbs[globalIdx]->lower.x = 1e6;
        aabbs[globalIdx]->lower.y = 1e6;
        aabbs[globalIdx]->lower.z = 1e6;

        aabbs[globalIdx]->upper.x = -1e6;
        aabbs[globalIdx]->upper.y = -1e6;
        aabbs[globalIdx]->upper.z = -1e6;

        int i1 = f[3*globalIdx + 0];
        int i2 = f[3*globalIdx + 1];
        int i3 = f[3*globalIdx + 2];


        aabbs[globalIdx]->lower.x = min(aabbs[globalIdx]->lower.x, start[3*i1 + 0]);
        aabbs[globalIdx]->upper.x = max(aabbs[globalIdx]->upper.x, start[3*i1 + 0]);
        aabbs[globalIdx]->lower.x = min(aabbs[globalIdx]->lower.x, end[3*i1 + 0]);
        aabbs[globalIdx]->upper.x = max(aabbs[globalIdx]->upper.x, end[3*i1 + 0]);

        aabbs[globalIdx]->lower.y = min(aabbs[globalIdx]->lower.y, start[3*i1 + 1]);
        aabbs[globalIdx]->upper.y = max(aabbs[globalIdx]->upper.y, start[3*i1 + 1]);
        aabbs[globalIdx]->lower.y = min(aabbs[globalIdx]->lower.y, end[3*i1 + 1]);
        aabbs[globalIdx]->upper.y = max(aabbs[globalIdx]->upper.y, end[3*i1 + 1]);

        aabbs[globalIdx]->lower.z = min(aabbs[globalIdx]->lower.z, start[3*i1 + 2]);
        aabbs[globalIdx]->upper.z = max(aabbs[globalIdx]->upper.z, start[3*i1 + 2]);
        aabbs[globalIdx]->lower.z = min(aabbs[globalIdx]->lower.z, end[3*i1 + 2]);
        aabbs[globalIdx]->upper.z = max(aabbs[globalIdx]->upper.z, end[3*i1 + 2]);

        aabbs[globalIdx]->lower.x = min(aabbs[globalIdx]->lower.x, start[3*i2 + 0]);
        aabbs[globalIdx]->upper.x = max(aabbs[globalIdx]->upper.x, start[3*i2 + 0]);
        aabbs[globalIdx]->lower.x = min(aabbs[globalIdx]->lower.x, end[3*i2 + 0]);
        aabbs[globalIdx]->upper.x = max(aabbs[globalIdx]->upper.x, end[3*i2 + 0]);

        aabbs[globalIdx]->lower.y = min(aabbs[globalIdx]->lower.y, start[3*i2 + 1]);
        aabbs[globalIdx]->upper.y = max(aabbs[globalIdx]->upper.y, start[3*i2 + 1]);
        aabbs[globalIdx]->lower.y = min(aabbs[globalIdx]->lower.y, end[3*i2 + 1]);
        aabbs[globalIdx]->upper.y = max(aabbs[globalIdx]->upper.y, end[3*i2 + 1]);

        aabbs[globalIdx]->lower.z = min(aabbs[globalIdx]->lower.z, start[3*i2 + 2]);
        aabbs[globalIdx]->upper.z = max(aabbs[globalIdx]->upper.z, start[3*i2 + 2]);
        aabbs[globalIdx]->lower.z = min(aabbs[globalIdx]->lower.z, end[3*i2 + 2]);
        aabbs[globalIdx]->upper.z = max(aabbs[globalIdx]->upper.z, end[3*i2 + 2]);

        aabbs[globalIdx]->lower.x = min(aabbs[globalIdx]->lower.x, start[3*i3 + 0]);
        aabbs[globalIdx]->upper.x = max(aabbs[globalIdx]->upper.x, start[3*i3 + 0]);
        aabbs[globalIdx]->lower.x = min(aabbs[globalIdx]->lower.x, end[3*i3 + 0]);
        aabbs[globalIdx]->upper.x = max(aabbs[globalIdx]->upper.x, end[3*i3 + 0]);

        aabbs[globalIdx]->lower.y = min(aabbs[globalIdx]->lower.y, start[3*i3 + 1]);
        aabbs[globalIdx]->upper.y = max(aabbs[globalIdx]->upper.y, start[3*i3 + 1]);
        aabbs[globalIdx]->lower.y = min(aabbs[globalIdx]->lower.y, end[3*i3 + 1]);
        aabbs[globalIdx]->upper.y = max(aabbs[globalIdx]->upper.y, end[3*i3 + 1]);

        aabbs[globalIdx]->lower.z = min(aabbs[globalIdx]->lower.z, start[3*i3 + 2]);
        aabbs[globalIdx]->upper.z = max(aabbs[globalIdx]->upper.z, start[3*i3 + 2]);
        aabbs[globalIdx]->lower.z = min(aabbs[globalIdx]->lower.z, end[3*i3 + 2]);
        aabbs[globalIdx]->upper.z = max(aabbs[globalIdx]->upper.z, end[3*i3 + 2]);

        printf("%i, %f , %f , %f , %f , %f , %f \n",globalIdx , aabbs[globalIdx]->upper.x,aabbs[globalIdx]->upper.y,aabbs[globalIdx]->upper.z,aabbs[globalIdx]->lower.x,aabbs[globalIdx]->lower.y,aabbs[globalIdx]->lower.z );

    }
}


// brute force broad phase
__global__ void simple_broad_phase(aabb** lhs, aabb* rhs, int* col, int ntris)
{


    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;


    //printf("RHS : %f , %f , %f , %f , %f , %f \n",rhs->lower.x , rhs->lower.y  , rhs->lower.z , rhs->upper.x , rhs->upper.y  , rhs->upper.z  );
    //printf("LHS : %f , %f , %f , %f , %f , %f \n",lhs[globalIdx].lower.x , lhs[globalIdx].lower.y  , lhs[globalIdx].lower.z , lhs[globalIdx].upper.x , lhs[globalIdx].upper.y  , lhs[globalIdx].upper.z  );

    if (globalIdx < ntris)
    {

        // printf("%f , %f , %f , %f , %f , %f \n",lhs[globalIdx].lower.x , lhs[globalIdx].lower.y  , lhs[globalIdx].lower.z , lhs[globalIdx].upper.x , lhs[globalIdx].upper.y  , lhs[globalIdx].upper.z  );

        if (lhs[globalIdx]->upper.x < rhs->lower.x || rhs->upper.x < lhs[globalIdx]->lower.x) {
            col[globalIdx] = false;
            //printf("%i - test 1 \n",globalIdx);
            return;
        }
        if (lhs[globalIdx]->upper.y < rhs->lower.y || rhs->upper.y < lhs[globalIdx]->lower.y) {
            col[globalIdx] = false;
            //printf("%i - test 2 \n",globalIdx);
            return;
        }
        if (lhs[globalIdx]->upper.z < rhs->lower.z || rhs->upper.z < lhs[globalIdx]->lower.z) {
            col[globalIdx] = false;
            //printf("%i - test 3 \n",globalIdx);
            return;
        }
    }

    col[globalIdx] = true;
    //printf("%i - collision\n",globalIdx);
    //printf("RHS : %f , %f , %f , %f , %f , %f \n",rhs->lower.x , rhs->lower.y  , rhs->lower.z , rhs->upper.x , rhs->upper.y  , rhs->upper.z  );
    //printf("LHS : %f , %f , %f , %f , %f , %f \n",lhs[globalIdx].lower.x , lhs[globalIdx].lower.y  , lhs[globalIdx].lower.z , lhs[globalIdx].upper.x , lhs[globalIdx].upper.y  , lhs[globalIdx].upper.z  );


}