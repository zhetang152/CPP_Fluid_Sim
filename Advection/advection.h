#pragma once
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\MACGrid.h"
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\grid.h"
#include "D:\Computation\FluidSim\CPP_Sim\Scene\Geometry\SolidShape.h"
namespace Advector {
    /**
     * @brief 半Lagrange平流方程
     * @param q_old 旧的标量场
     * @param velocitygrid 速度网格
     * @param dt 时间步长
     * @return 新的标量场 
     */
    Grid<float> advect(const Grid<float>& q_old, const MACGrid& velocityGrid, float dt);
    void advect_velocity(MACGrid& grid, float dt);
    void advect_particles(MACGrid& grid, const SolidShape& solid, float dt);
}