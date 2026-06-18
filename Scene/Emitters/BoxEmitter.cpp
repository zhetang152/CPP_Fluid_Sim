#include "BoxEmitter.h"

#include <cstdlib>

void BoxEmitter::emit(MACGrid& grid, float dt) {
    float num_new_particles_float = m_particles_per_second * dt;
    int num_new_particles_int = static_cast<int>(num_new_particles_float);
    float fractional_part = num_new_particles_float - num_new_particles_int;

    if (static_cast<float>(rand()) / RAND_MAX < fractional_part) {
        ++num_new_particles_int;
    }

    for (int i = 0; i < num_new_particles_int; ++i) {
        float rand_x = static_cast<float>(rand()) / RAND_MAX;
        float rand_y = static_cast<float>(rand()) / RAND_MAX;
        float rand_z = static_cast<float>(rand()) / RAND_MAX;

        Point3f particle_pos = m_min_corner + Vector3f(
            rand_x * (m_max_corner.x - m_min_corner.x),
            rand_y * (m_max_corner.y - m_min_corner.y),
            rand_z * (m_max_corner.z - m_min_corner.z));

        Particles new_particle;
        new_particle.position = particle_pos;
        new_particle.velocity = m_initial_velocity;
        grid.particles().push_back(new_particle);
    }
}
