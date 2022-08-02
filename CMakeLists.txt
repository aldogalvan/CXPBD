cmake_minimum_required(VERSION 3.21)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_COMPILER /usr/local/cuda/bin/nvcc)
set(CMAKE_CUDA_STANDARD 14)

project(CXPBD CXX CUDA)

#CHAI-3D
# include Chai3d
# -- this is required by Chai3d to pick the correct path for OpenGL headers

FIND_PACKAGE(CHAI3D REQUIRED)
ADD_DEFINITIONS(${CHAI3D_DEFINITIONS})
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
    ADD_DEFINITIONS(-D_MACOSX -DMACOSX)
elseif (CMAKE_SYSTEM_NAME MATCHES Linux)
    ADD_DEFINITIONS(-D_LINUX -DLINUX)
endif () #NOTE: Windows is currently unsupported
LINK_DIRECTORIES(${CHAI3D_LIBRARY_DIRS})

# Add the include directory to the include paths
INCLUDE_DIRECTORIES (
        ${CHAI3D_INCLUDE_DIRS}
)

#

INCLUDE_DIRECTORIES (
        src
        src/world
        src/collision
        src/constraints
)

# static library source files
file (GLOB_RECURSE CXPBD_SOURCES src/*.cpp src/*.h)

# static library
add_library (CXPBD STATIC ${CXPBD_SOURCES} src/collision/CXPBDContactConstraint.cpp src/collision/CXPBDContactConstraint.h src/collision/CXPBDContactProjection.cpp src/collision/CXPBDContactProjection.h src/collision/CXPBDProxyAlgorithm.cpp src/collision/CXPBDProxyAlgorithm.h src/collision/CXPBDEdgeCollisionDetection.cpp src/collision/CXPBDEdgeCollisionDetection.h)

# use cuda
if (CXPBD_USE_CUDA)
    find_package(CUDA REQUIRED)
    set(CMAKE_CUDA_ARCH 5 6)
    set(CMAKE_CONFIGURATION_TYPES "Release")
endif()



#include LibIGL
add_subdirectory(${PROJECT_SOURCE_DIR}/external/libigl-main)
include_directories(${PROJECT_SOURCE_DIR}/external/libigl-main)

# GLFW
if (EXISTS ${PROJECT_SOURCE_DIR}/extras/GLFW)
    add_subdirectory (${PROJECT_SOURCE_DIR}/extras/GLFW)
endif ()

include_directories(extras/GLFW/include)

#EIGEN3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(CXPBD
        PUBLIC
        ${CHAI3D_LIBRARIES}
        ${GLFW_LIBRARIES}
        Eigen3::Eigen
        igl::core
        )

add_subdirectory(examples)
