cmake_minimum_required(VERSION 3.21)
project(CXPBD)

set(CMAKE_CXX_STANDARD 14)


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

INCLUDE_DIRECTORIES (
        src
)


# static library source files
file (GLOB_RECURSE CXPBD_SOURCES src/*.cpp src/*.h)

# static library
add_library (CXPBD STATIC ${CXPBD_SOURCES} src/CXPBDTool.cpp src/CXPBDTool.h src/CXPBDAABB.cpp src/CXPBDAABB.h)

#include LibIGL
target_include_directories(CXPBD PUBLIC ${PROJECT_SOURCE_DIR}/external/libigl-main )

# GLFW
if (EXISTS ${PROJECT_SOURCE_DIR}/extras/GLFW)
    add_subdirectory (${PROJECT_SOURCE_DIR}/extras/GLFW)
endif ()

#EIGEN3
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

target_link_libraries(CXPBD
        PUBLIC
        ${CHAI3D_LIBRARIES}
        ${GLFW_LIBRARIES}
        Eigen3::Eigen
        )

add_subdirectory(examples)
