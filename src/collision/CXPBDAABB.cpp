//
// Created by aldof on 4/13/2022.
//

#include <iostream>
#include "CXPBDAABB.h"

using namespace std;

struct BBox{
    double mins[3];
    double maxs[3];
};

bool intersects(const BBox& b1, const BBox& b2)
{
    for (int i = 0; i < 3; i++) {
        if ((b1.maxs[i] < b2.mins[i]) || (b2.maxs[i] < b1.mins[i]))
            return false;
    }
    return true;
}

struct AABBNode {
    AABBNode()
    {
    }
    ~AABBNode()
    {
    }

    std::shared_ptr<AABBNode> left;
    std::shared_ptr<AABBNode> right;
    BBox box;
    int childtri = -1;

};

class NodeComparator {

public:
    NodeComparator(int axis)
            : axis(axis)
    {
    }

    int axis;

    bool operator()(const std::shared_ptr<const AABBNode>& left,
                    const std::shared_ptr<const AABBNode>& right) const
    {
        return left->box.mins[axis] < right->box.mins[axis];
    }
};

std::shared_ptr<AABBNode>
buildAABB(vector<std::shared_ptr<AABBNode>>& nodes)
{

    if (nodes.size() == 0)
        return NULL;
    else if (nodes.size() == 1)
        return nodes[0];

    double axismins[3];
    double axismaxs[3];
    for (int i = 0; i < 3; i++) {
        axismins[i] = numeric_limits<double>::infinity();
        axismaxs[i] = -numeric_limits<double>::infinity();
    }
    int nnodes = (int)nodes.size();
    for (int i = 0; i < nnodes; i++) {
        for (int j = 0; j < 3; j++) {
            axismins[j] = min(axismins[j], nodes[i]->box.mins[j]);
            axismaxs[j] = max(axismaxs[j], nodes[i]->box.maxs[j]);
        }
    }
    double widths[3];
    for (int i = 0; i < 3; i++)
        widths[i] = axismaxs[i] - axismins[i];
    int splitaxis = -1;
    if (widths[0] >= widths[1] && widths[0] >= widths[2])
        splitaxis = 0;
    else if (widths[1] >= widths[0] && widths[1] >= widths[2])
        splitaxis = 1;
    else
        splitaxis = 2;
    std::sort(nodes.begin(), nodes.end(), NodeComparator(splitaxis));
    vector<std::shared_ptr<AABBNode>> left(nnodes / 2);
    vector<std::shared_ptr<AABBNode>> right(nnodes - nnodes / 2);
    for (int i = 0; i < nnodes / 2; i++) {
        left[i] = nodes[i];
    }
    for (int i = nnodes / 2; i < nnodes; i++) {
        right[i - nnodes / 2] = nodes[i];
    }
    auto node = std::make_shared<AABBNode>();
    node->left = buildAABB(left);
    node->right = buildAABB(right);
    for (int i = 0; i < 3; i++) {
        node->box.mins[i] = min(node->left->box.mins[i], node->right->box.mins[i]);
        node->box.maxs[i] = max(node->left->box.maxs[i], node->right->box.maxs[i]);
    }
    return node;
}

std::shared_ptr<AABBNode> buildAABB(Eigen::Ref<const Eigen::MatrixXd> startPos,
                                    Eigen::Ref<const Eigen::MatrixXd> endPos,
                                    Eigen::Ref<const Eigen::MatrixXi> F)
{
    int ntris = F.rows();
    vector<std::shared_ptr<AABBNode>> leaves(ntris);
    for (int j = 0; j < ntris; j++) {
        auto leaf = std::make_shared<AABBNode>();
        leaf->childtri = j;
        Eigen::Vector3i face = F.row(j).transpose();
        BBox box;
        for (int k = 0; k < 3; k++) {
            box.mins[k] = numeric_limits<double>::infinity();
            box.maxs[k] = -numeric_limits<double>::infinity();
        }
        for (int k = 0; k < 3; k++) {
            Eigen::Vector3d startpoint = startPos.row(face[k]).transpose();
            Eigen::Vector3d endpoint = endPos.row(face[k]).transpose();
            for (int l = 0; l < 3; l++) {
                box.mins[l] = min(box.mins[l], startpoint[l]);
                box.maxs[l] = max(box.maxs[l], startpoint[l]);
                box.mins[l] = min(box.mins[l], endpoint[l]);
                box.maxs[l] = max(box.maxs[l], endpoint[l]);
            }
        }
        leaf->box = box;
        leaves[j] = leaf;
    }
    return buildAABB(leaves);
}


std::shared_ptr<AABBNode> buildAABB(Eigen::Ref<const Eigen::MatrixXd> pos, Eigen::Ref<const Eigen::MatrixXi> F )
{
    int ntris = F.rows();
    vector<std::shared_ptr<AABBNode>> leaves(ntris);
    for (int j = 0; j < ntris; j++) {
        auto leaf = std::make_shared<AABBNode>();
        leaf->childtri = j;
        Eigen::Vector3i face = F.row(j).transpose();
        BBox box;
        for (int k = 0; k < 3; k++) {
            box.mins[k] = numeric_limits<double>::infinity();
            box.maxs[k] = -numeric_limits<double>::infinity();
        }
        for (int k = 0; k < 3; k++) {
            Eigen::Vector3d point = pos.row(face[k]).transpose();
            for (int l = 0; l < 3; l++) {
                box.mins[l] = min(box.mins[l], point[l]);
                box.maxs[l] = max(box.maxs[l], point[l]);
            }
        }
        leaf->box = box;
        leaves[j] = leaf;
    }
    return buildAABB(leaves);
}


std::shared_ptr<BBox> buildAABB(Eigen::Ref<const Eigen::Vector3d> startPos,
                                    Eigen::Ref<const Eigen::Vector3d> endPos,
                                    double radius)
{

    auto box = make_shared<BBox>();

    for (int k = 0; k < 3; k++)
    {
        box->mins[k] = numeric_limits<double>::infinity();
        box->maxs[k] = -numeric_limits<double>::infinity();
    }

    for (int l = 0; l < 3; l++) {
        box->mins[l] = min(box->mins[l], startPos(l));
        box->maxs[l] = max(box->maxs[l], startPos(l));
        box->mins[l] = min(box->mins[l], endPos(l)) - radius;
        box->maxs[l] = max(box->maxs[l], endPos(l)) + radius;
    }


    return box;
}


std::shared_ptr<BBox> buildAABB(Eigen::Ref<const Eigen::Vector3d> pos,
                                double radius)
{
    auto box = make_shared<BBox>();

    for (int k = 0; k < 3; k++)
    {
        box->mins[k] = numeric_limits<double>::infinity();
        box->maxs[k] = -numeric_limits<double>::infinity();
    }

    for (int l = 0; l < 3; l++) {
        box->mins[l] = min(box->mins[l], pos(l)) - radius;
        box->maxs[l] = max(box->maxs[l], pos(l)) + radius;
    }


    return box;
}

void intersect(std::shared_ptr<AABBNode> left,
               std::shared_ptr<AABBNode> right,
               std::vector<Collision>& collisions)
{

    if (!intersects(left->box, right->box))
        return;
    if (left->childtri != -1)
    {
        if (right->childtri != -1)
        {
            collisions.emplace_back(Collision{left->childtri, right->childtri});
        } else {
            intersect(left, right->left, collisions);
            intersect(left, right->right, collisions);
        }
    } else
    {
        intersect(left->left, right, collisions);
        intersect(left->right, right, collisions);
    }
}

void intersect(std::shared_ptr<BBox> left,
               std::shared_ptr<AABBNode> right,
               std::vector<Collision>& collisions)
{


    if (!intersects(*left, right->box))
        return;

    if (right->childtri != -1)
    {
        collisions.emplace_back(Collision{NULL, right->childtri});
    } else {
        intersect(left, right->left, collisions);
        intersect(left, right->right, collisions);
    }
}