//
// Created by aldof on 4/13/2022.
//

#ifndef CXPBD_CXPBDAABB_H
#define CXPBD_CXPBDAABB_H

#include <Eigen/Core>
#include <memory>
#include <vector>



struct AABBNode;
struct BBox;

std::shared_ptr<AABBNode> buildAABB(Eigen::Ref<const Eigen::MatrixXd> startPos,
                                    Eigen::Ref<const Eigen::MatrixXd> endPos,
                                    Eigen::Ref<const Eigen::MatrixXi> F);

std::shared_ptr<BBox> buildAABB(Eigen::Ref<const Eigen::Vector3d> startPos,
                        Eigen::Ref<const Eigen::Vector3d> endPos , double radius);

struct Collision {
    int collidingTriangle1; // Triangle from the "left" AABBNode
    int collidingTriangle2; // Triangle from the "right" AABBNode
};

// Intersect AABB "left" vs "right"
void intersect(std::shared_ptr<AABBNode> left,
               std::shared_ptr<AABBNode> right,
               std::vector<Collision>& collisions);

// Intersect BB and AABB
void intersect(std::shared_ptr<BBox> left,
               std::shared_ptr<AABBNode> right,
               std::vector<Collision>& collisions);




#endif //CXPBD_CXPBDAABB_H
