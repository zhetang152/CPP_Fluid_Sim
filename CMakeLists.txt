cmake_minimum_required(VERSION 3.10)
project(CPP_SIM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 添加项目的根目录到头文件搜索路径
include_directories(${CMAKE_SOURCE_DIR})

# 源文件列表
set(SOURCES
    main.cpp

    Grid_Construction/vector3D.cpp
    Grid_Construction/SDFUtils.cpp

    Boundary/SolidBoundary.cpp

    Scene/SceneManager.cpp
    Scene/Emitters/BoxEmitter.cpp
    Scene/Surface/MarchingCubes.cpp

    Solver/solver.cpp

    Advection/advection.cpp
    Advection/FLIP.cpp

    IO/DataExporter.cpp
    
    Force/General_Force/GravityForce.cpp
)
# 生成可执行文件
add_executable(main ${SOURCES})