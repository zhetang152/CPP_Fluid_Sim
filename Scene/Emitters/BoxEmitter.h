#pragma once

#include "Grid_Construction/vecmath.h"
#include "ParticleEmitter.h"

class BoxEmitter : public ParticleEmitter {
private:
    Point3f m_min_corner;
    Point3f m_max_corner;
    float m_particles_per_second;
    Vector3f m_initial_velocity;

public:
    BoxEmitter(
        const Point3f& min_corner,
        const Point3f& max_corner,
        float particles_per_second,
        const Vector3f& initial_velocity)
        : m_min_corner(min_corner),
          m_max_corner(max_corner),
          m_particles_per_second(particles_per_second),
          m_initial_velocity(initial_velocity) {}

    void emit(MACGrid& grid, float dt) override;
};
