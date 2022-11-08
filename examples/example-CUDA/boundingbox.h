#ifndef CUDAISFUN_BOUNDINGBOX_H
#define CUDAISFUN_BOUNDINGBOX_H

#include <vector>

struct BBox
{
    float mins[3];
    float maxs[3];
};

struct AABBNode
{
    AABBNode* left;
    AABBNode* right;
    BBox box;
    int childtri = -1;
};

// this function will build a full bounding box for an object mesh
AABBNode* buildAABB(const float* startPos, const float* endPos, const int* F, const int ntris);

// this function will build a bounding box for a point primitive
BBox* buildBB(const float* startPos, const float* endPos);

// this function computes an intersection point
void intersect(const BBox* left, const AABBNode* right, std::vector<int>& collisions);

#endif //CUDAISFUN_BOUNDINGBOX_H