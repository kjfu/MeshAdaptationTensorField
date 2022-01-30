cmake_minimum_required(VERSION 3.5)
project(MeshAdaptationTensorField)
FIND_PACKAGE( OpenMP REQUIRED)

if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(
    ./src
    ./external/eigen
)
set(MESH_ADAPTATION_TENSOR_FIELD_SRCS
    src/main.cpp
    src/mesh.cpp
    src/parser.cpp
)

add_executable(MeshAdaptationTensorField ${MESH_ADAPTATION_TENSOR_FIELD_SRCS})