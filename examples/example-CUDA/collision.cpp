#include "collision.h"

// TODO: Best way to get collision info
void findCollisions(const float* d_goal, float* d_proxy, meshObject* obj)
{

    // broad phase collision
    auto potCol = CTCD::broadPhase(d_goal, d_proxy, obj);


    // narrow phase collision detection
    CTCD::narrowPhase(d_goal,d_proxy,obj,potCol);

}

vector<int> CTCD::broadPhase(const float* d_goal, float* d_proxy, meshObject* obj)
{

    // computes the normals for the object
    obj->computeNormals();

    // Builds a bounding box for the tool
    auto toolBB_ = buildBB(d_goal,d_proxy);

    // the vector containing the possible collisions
    std::vector<int> potentialCollisions;

    // Broad phase collision detection
    intersect(toolBB_, obj->buildBoundingBox(), potentialCollisions);

    return potentialCollisions;
}

