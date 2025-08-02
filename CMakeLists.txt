cmake_minimum_required(VERSION 3.10)
project(CPP_SIM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 源文件列表
include_directories(${CMAKE_SOURCE_DIR})
set(SOURCES
    main.cpp
    Grid_Construction/vector3D.cpp
    Boundary/SolidBoundary.cpp
    Scene/SceneManager.cpp
    Solver/solver.cpp
)

add_executable(main ${SOURCES})