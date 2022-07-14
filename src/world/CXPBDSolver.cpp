//
// Created by aldof on 4/2/2022.
//


#include "CXPBDSolver.h"


void solve(
        cXPBDDeformableMesh* model,
        cXPBDToolMesh* tool,
        Eigen::MatrixX3d const& fext,
        double timestep,
        std::uint32_t iterations,
        std::uint32_t substeps,
        bool gravityEnabled
        )
{

    auto const num_iterations = iterations / substeps;
    double dt                 = timestep / static_cast<double>(substeps);
    auto const& constraints   = model->constraints();
    auto const J              = constraints.size();
    std::vector<double> lagrange_multipliers(J, 0.);

    for (auto s = 0u; s < substeps; ++s)
    {
        auto& v = model->velocity();
        auto& x = model->positions();

        tool->updatePos();

        // step the graphics forward (test????)
        // tool->setPnext(dt);
        // tool->buildAABBBoundaryBox();

        auto const& m = model->mass();
        Eigen::MatrixX3d a = fext.array().colwise() / m.array();

        if (gravityEnabled)
        {
            Eigen::RowVector3d g(0,0,-9.8);
            a.rowwise() += g;
        }

        // explicit euler step
        auto vexplicit = v + dt * a;
        Eigen::MatrixXd p = x + dt * vexplicit;

        // generate collision constraints here ...
        double t;
        bool collision = false;
        Eigen::Vector3d pos;
        double t_;
        std::vector<ColInfo*> collisions;
        collision = findCollisions(pos,pos,0.1,model,collisions);

        if (collision == true)
        {
            //applyVertexVsBodyConstraints(collisions, p, model, tool);
            tool->projectPos(t);
        }

        for (auto i = 0u ; i < model->numVerts() ; i++)
        {
            if (p(i,2) < -.1)
            {
                p(i,2) = -.1;
            }
        }

        Eigen::Vector3d F(0,0,0);
        // sequential gauss seidel type solve
        std::fill(lagrange_multipliers.begin(), lagrange_multipliers.end(), 0.0);
        for (auto n = 0u; n < num_iterations; ++n)
        {
            for (auto j = 0u; j < J; ++j)
            {
                auto const& constraint = constraints[j];
                constraint->project(p, x, m, lagrange_multipliers[j], dt, F);
            }
        }

        double lagrange_sum = 0;
        
        for (auto it = lagrange_multipliers.begin() ; it != lagrange_multipliers.end(); it++)
        {
            lagrange_sum += *it;
        }

        // update solution
        for (auto i = 0u; i < x.rows(); ++i)
        {
            v.row(i) = (p.row(i) - x.row(i)) / dt;
            x.row(i) = p.row(i);
        }

        model->updateChai3d();
        tool->updateChai3d();
    }
}

