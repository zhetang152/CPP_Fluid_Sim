#pragma once
#include "Grid_Construction\GridAndParticleSystem.h"
#include "Grid_Construction\grid.h"
#include "Scene\Geometry\SolidShape.h"
namespace Advector {
    /**
     * @brief 半Lagrange平流方程
     * @param q_old 旧的标量场
     * @param velocitygrid 速度网格
     * @param dt 时间步长
     * @return 新的标量场 
     */
    //底层重载，直接接收网格引用，避免对象构造开销
    Vector3f get_velocity_at(const Grid<float>& u, const Grid<float>& v, const Grid<float>& w, float dx, const Point3f& pos);
    //便捷接口，内部调用上面的重载
    Vector3f get_velocity_at(const MACGrid& grid, const Point3f& pos);
    Grid<float> advect(const Grid<float>& q_old, const MACGrid& velocityGrid, float dt);
    void advect_velocity(MACGrid& grid, float dt);
    void advect_particles(MACGrid& grid, const SolidShape& solid, float dt);
    
}