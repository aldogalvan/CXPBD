foreach (example example-01)
    file (GLOB source ${example}/*.cpp ${example}/*.h)
    #set(tetgen_DIR "/home/aldo/tetgen1.6.0/build")
    #find_package(tetgen)
    include_directories(../src)
    #include_directories(/home/aldo/tetgen-1.6.0/build)
    add_executable (${example} ${source})
    target_link_libraries (${example} CXPBD )
    add_definitions (-DTETLIBRARY)
endforeach ()


#foreach (example example-CUDA)
 #   find_package(CUDA)
  #  INCLUDE_DIRECTORIES(
   #         ${CUDA_SDK_ROOT_DIR}/C/common/inc/
    #)
    #file (GLOB source ${example}/*.cu ${example}/*.cpp ${example}/*.h)
    #include_directories(../src)
    #add_executable(${example} ${source} example-CUDA/example-CUDA.cpp)
    #target_link_libraries (${example} CXPBD cuda)
    #add_definitions (-DTETLIBRARY)
#endforeach()
