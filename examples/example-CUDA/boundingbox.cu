#include "boundingbox_kernel.cuh"
#include "boundingbox.h"
#include "thrust/device_vector.h"
#include "thrust/host_vector.h"

using namespace std;




struct NodeComparator {

    __host__ __device__ NodeComparator(int axis) : axis(axis)
    {
    }

    int axis;

    __host__ __device__
    bool operator()(const AABBNode* left,
                    const AABBNode* right) const
    {
        //printf("%f \n",left->box.mins[axis]);
        //printf("HEY \n");
        return left->box.mins[axis] < right->box.mins[axis];
    }
};


void intersects(const BBox& b1, const BBox& b2, bool col)
{
    for (int i = 0; i < 3; i++)
    {
        if ((b1.maxs[i] < b2.mins[i]) || (b2.maxs[i] < b1.mins[i]))
            col = 0;
    }
    col = 1;
}

AABBNode* buildAABB(thrust::host_vector<AABBNode*> nodes)
{

    if (nodes.size() == 0)
        return NULL;
    else if (nodes.size() == 1)
        return nodes[0];

    float axismins[3];
    float axismaxs[3];

    for (int i = 0; i < 3; i++)
    {
        axismins[i] = INFINITY;
        axismaxs[i] = -INFINITY;
    }

    int nnodes = nodes.size();

    for (int i = 0; i < nnodes; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            axismins[j] = min(axismins[j], nodes[i]->box.mins[j]);
            axismaxs[j] = max(axismaxs[j], nodes[i]->box.maxs[j]);
        }
    }

    //thrust::device_vector<AABBNode*> d_nodes = nodes;

    float widths[3];
    for (int i = 0; i < 3; i++)
        widths[i] = axismaxs[i] - axismins[i];
    int splitaxis = -1;

    if (widths[0] >= widths[1] && widths[0] >= widths[2])
        splitaxis = 0;
    else if (widths[1] >= widths[0] && widths[1] >= widths[2])
        splitaxis = 1;
    else
        splitaxis = 2;

    // TODO: make sure to check whether raw pointer works better
    //std::cout << "ONE" << std::endl;
    // thrust::sort( thrust::device,d_nodes.begin(), d_nodes.end(),NodeComparator(splitaxis));
    //cudaDeviceSynchronize();
    //std::cout << "TWO" << std::endl;
    thrust::sort(nodes.begin(),nodes.end(), NodeComparator(splitaxis));

    int nnodesleft =  nnodes / 2;
    int nnodesright = nnodes - nnodes / 2;

    // nodes = d_nodes;
    thrust::host_vector<AABBNode*> left(nnodesleft);
    thrust::host_vector<AABBNode*> right(nnodesright);

    thrust::copy(nodes.begin(),nodes.begin() + nnodesleft, left.begin());
    thrust::copy(nodes.begin() + nnodesleft ,nodes.end(), right.begin());

    AABBNode* node = new AABBNode;


    node->left = buildAABB(left);
    node->right = buildAABB(right);

    for (int i = 0; i < 3; i++)
    {
        node->box.mins[i] = min(node->left->box.mins[i], node->right->box.mins[i]);
        node->box.maxs[i] = max(node->left->box.maxs[i], node->right->box.maxs[i]);
    }

    return node;
}

AABBNode* buildAABB(const float* startPos, const float* endPos, const int* F, const int ntris)
{

    thrust::device_vector<AABBNode*> leaves(ntris);
    AABBNode** pleaves = thrust::raw_pointer_cast(leaves.data());
    buildAABB<<<2,1024>>>(startPos,endPos,F,pleaves,ntris);

    thrust::host_vector<AABBNode*> h_leaves = leaves;

    return buildAABB(h_leaves);
}

BBox* buildBB(const float* startPos, const float* endPos)
{
    BBox* bbox = new BBox;

    for (int k = 0; k < 3; k++)
    {
        bbox->mins[k] = INFINITY;
        bbox->maxs[k] = -INFINITY;
    }

    for (int l = 0; l < 3; l++) {
        bbox->mins[l] = min(bbox->mins[l], startPos[l]);
        bbox->maxs[l] = max(bbox->maxs[l], startPos[l]);
        bbox->mins[l] = min(bbox->mins[l], endPos[l]);
        bbox->maxs[l] = max(bbox->maxs[l], endPos[l]);
    }

    return bbox;
}

void intersect(const BBox* left, const AABBNode* right, vector<int>& collisions)
{

    bool pBoolean;
    intersects(*left, right->box,pBoolean);

    if (right->childtri != -1)
    {
        collisions.push_back(right->childtri);
    } else {
        intersects(*left, right->box,pBoolean);
        intersects(*left, right->box,pBoolean);
    }
}
