#pragma once
#include "Grid_Construction/vector3D.hpp"
#include "ParticleEmitter.hpp"

class BoxEmitter : public ParticleEmitter{
private:
    Vector3D m_min_corner;// 长方体的一个角点
    Vector3D m_max_corner;// 长方体的对角点
    float m_particles_per_second;//每秒发射的粒子数
    Vector3D m_initial_velocity;//发射粒子的初速度
public:
    BoxEmitter(
        const Vector3D& min_corner,
        const Vector3D& max_corner,
        float particles_per_second,
        const Vector3D& initial_velocity
    ):m_min_corner(min_corner), m_max_corner(max_corner), 
      m_particles_per_second(particles_per_second), m_initial_velocity(initial_velocity) {}
    void emit(MACGrid& grid, float dt) override;
};