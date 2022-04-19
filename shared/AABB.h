#ifndef PSIM_CORE_CLOTH2_AABB_H
#define PSIM_CORE_CLOTH2_AABB_H

#include <Eigen/Core>
#include <memory>
#include <vector>

struct AABBNode;

std::shared_ptr<AABBNode> buildAABB(Eigen::Ref<const Eigen::MatrixXd> startPos,
                                    Eigen::Ref<const Eigen::MatrixXd> endPos,
                                    Eigen::Ref<const Eigen::MatrixXi> F);

struct Collision {
    int collidingTriangle1; // Triangle from the "left" AABBNode
    int collidingTriangle2; // Triangle from the "right" AABBNode
};

// Intersect AABB "left" vs "right"
// CAVEAT: this function will *NOT* clear the content of "collisions" for you.
void intersect(std::shared_ptr<AABBNode> left,
               std::shared_ptr<AABBNode> right,
               std::vector<Collision>& collisions);



#endif
