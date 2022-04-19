#include "helper.h"
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>

std::tuple<Eigen::MatrixX3d, Eigen::MatrixX3i>
loadOBJ(const std::string& fn)
{
        Eigen::MatrixX3d V;
        Eigen::MatrixX3i F;
        igl::readOBJ(fn, V, F);
        return std::make_tuple(V, F);
}

void saveOBJ(const Eigen::Matrix<double, -1, -1>& V,
             const Eigen::Matrix<int, -1, -1>& F,
             const std::string& fn)
{
    igl::writeOBJ(fn, V, F);
}

