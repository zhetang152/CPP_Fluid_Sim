cmake_minimum_required(VERSION 3.10)
project(CPP_SIM)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(${CMAKE_SOURCE_DIR})

set(SOURCES
    main.cpp

    Grid_Construction/SDFUtils.cpp

    Boundary/SolidBoundary.cpp

    Scene/SceneManager.cpp
    Scene/Emitters/BoxEmitter.cpp
    Scene/Surface/ply.c
    Scene/Surface/MarchingCubes.cpp

    Solver/solver.cpp

    Advection/advection.cpp
    Advection/FLIP.cpp

    IO/DataExporter.cpp

    Force/General_Force/GravityForce.cpp
)

add_executable(main ${SOURCES})

if(MSVC)
    target_compile_options(main PRIVATE /utf-8)
endif()
