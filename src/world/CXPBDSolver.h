//
// Created by aldof on 4/2/2022.
//

#ifndef CXPBD_CXPBDSOLVER_H
#define CXPBD_CXPBDSOLVER_H

#include <Eigen/Dense>
#include "CXPBDDeformableObject.h"
#include "collision/CXPBDContinuousCollisionDetection.h"
#include "CXPBDSolver.h"
#include "CXPBDToolMesh.h"
#include "CXPBDTool.h"
#include "chai3d.h"

using namespace chai3d;

    void solve(
            cXPBDDeformableMesh* model,
            cXPBDToolMesh* tool,
            Eigen::MatrixX3d const& fext,
            double dt                = 0.01,
            std::uint32_t iterations = 10,
            std::uint32_t substeps   = 10,
            bool gravityEnabled = false);




#endif //CXPBD_CXPBDSOLVER_H
