#include "boundingbox.h"
#include <Eigen/Dense>

using namespace Eigen;

__global__ void buildAABB(const float* startPos, const float* endPos, const int* F, AABBNode** leaves, int ntris)
{


    int globalIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalIdx < ntris)
    {
        AABBNode* leaf = new AABBNode;
        leaf->childtri = globalIdx;
        Vector3i face(F[3 * globalIdx], F[3 * globalIdx + 1], F[3 * globalIdx + 2]);
        BBox box;

        for (int k = 0; k < 3; k++) {
            box.mins[k] = INFINITY;
            box.maxs[k] = -INFINITY;
        }

        for (int k = 0; k < 3; k++) {
            Vector3f startpoint(startPos[3 * face(k)], startPos[3 * face(k) + 1], startPos[3 * face(k) + 2]);
            Vector3f endpoint(endPos[3 * face(k)], endPos[3 * face(k) + 1], endPos[3 * face(k) + 2]);

            for (int l = 0; l < 3; l++) {
                box.mins[l] = min(box.mins[l], startpoint[l]);
                box.maxs[l] = max(box.maxs[l], startpoint[l]);
                box.mins[l] = min(box.mins[l], endpoint[l]);
                box.maxs[l] = max(box.maxs[l], endpoint[l]);

            }
        }

        leaf->box = box;
        leaves[globalIdx] = leaf;

    }

}

__global__ void aabbMinMax(AABBNode** nodes, float* axismins, float* axismaxs, int nnodes)
{
    for (int i = 0; i < 3; i++)
    {
        axismins[i] = INFINITY;
        axismaxs[i] = -INFINITY;
    }

    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            axismins[j] = min(axismins[j], nodes[i]->box.mins[j]);
            axismaxs[j] = max(axismaxs[j], nodes[i]->box.maxs[j]);

        }
    }

}