foreach (example example-01 )
    file (GLOB source ${example}/*.cpp ${example}/*.h)
    include_directories(../src)
    add_executable (${example} ${source})
    target_link_libraries (${example} CXPBD )
    add_definitions (-DTETLIBRARY)
endforeach ()

foreach (example example-CUDA)
    find_package(CUDA)
    INCLUDE_DIRECTORIES(
            ${CUDA_SDK_ROOT_DIR}/C/common/inc/
    )

    file (GLOB source ${example}/*.cu ${example}/*.cuh ${example}/*.cpp ${example}/*.h ${example}/lbvh/*.cuh
            ${example}/tetgen/*.h ${example}/tetgen/*.cpp)
    add_executable(${example} ${source})
    target_link_libraries (${example} CXPBD)
    add_subdirectory(${example}/Tight-Inclusion)
    add_definitions (-DTETLIBRARY)
endforeach()
