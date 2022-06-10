//
// Created by aldof on 4/2/2022.
//


#include "CXPBDSolver.h"

/*
void applyVertexVsBodyConstraints(const ColInfo& col_info, Eigen::MatrixXd& p,
                                  cXPBDDeformableMesh* model,
                                  cXPBDToolMesh* tool)
{
    const auto& vertCollisions = col_info.vertCollisions;
    const auto& p_tool = tool->positions();
    const auto& p_tool_last = tool->positions_last();
    // const auto& p_tool = tool->positions_next();
    const auto& f_def = model->faces();
    // const auto& f_tool = tool->faces();

    for (int i = 0; i < tool->numVerts(); i++){
        if (!vertCollisions[i].empty()){
            for (auto &col : vertCollisions[i]){
                Eigen::Vector3d pos = p_tool.row(i) ;
                Eigen::Vector3i face = f_def.row(col);
                Eigen::Vector3d p1 = p.row(face(0));
                Eigen::Vector3d p2 = p.row(face(1));
                Eigen::Vector3d p3 = p.row(face(2));
                auto p0 = (p1 + p2 + p3) / 3;
                Eigen::Vector3d norm = (p2 - p1).cross(p3 - p1).normalized();
                double j = (pos - p0).dot(norm);

                // std::cout << j << "," << norm.transpose() << std::endl;

                if ( j < 0 )
                {
                    p.row(face(0)) -=  2 * j * norm.transpose();
                    p.row(face(1)) -=  2 * j * norm.transpose();
                    p.row(face(2)) -=  2 * j * norm.transpose();
                }
                else
                {
                    p.row(face(0)) +=  2 * j * norm.transpose();
                    p.row(face(1)) +=  2 * j * norm.transpose();
                    p.row(face(2)) +=  2 * j * norm.transpose();
                }

                // std::cout << j << std::endl;
                //p.row(face(0)) +=  2 * j * norm.transpose();
                //p.row(face(1)) +=  2 * j * norm.transpose();
                //p.row(face(2)) +=  2 * j * norm.transpose();


            }
        }
    }
}

void applyFaceVsVertexConstraints(const ColInfo& col_info, Eigen::MatrixXd& p,
                                  cXPBDDeformableMesh* model,
                                  cXPBDToolMesh* tool)
{

    const auto& faceCollisions = col_info.faceCollisions;
    auto& p_tool = tool->positions();
    // const auto& p_tool = tool->positions_next();
    auto& f_def = model->faces();
    auto& f_tool = tool->faces();
    for (int i = 0; i < tool->numFaces(); i++){
        if (!faceCollisions[i].empty()){
            for (auto &col : faceCollisions[i]){
                Eigen::Vector3i face = f_tool.row(i);
                Eigen::Vector3d p1 = p_tool.row(face(0));
                Eigen::Vector3d p2 = p_tool.row(face(1));
                Eigen::Vector3d p3 = p_tool.row(face(2));
                auto p0 = (p1+p2+p3)/3;
                Eigen::Vector3d pos = p.row(col);
                Eigen::Vector3d norm = (p2 - p1).cross(p3 - p1).normalized();
                double j = (pos - p0).dot(norm);
                if (j < 0) {
                    p.row(i) -=  j * norm.transpose();
                }
                else{
                    //p.row(i) -=  j * norm.transpose();
                }

            }
        }
    }
}

*/

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
        const auto& collisions = findCollisions(p,pos,pos,*model,t_);

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

        model->updateChai3d(p);
        tool->updateChai3d();
    }
}

