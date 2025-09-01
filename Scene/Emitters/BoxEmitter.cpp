#include "BoxEmitter.hpp"
#include <cstdlib>//rand()
#include <vector>

void BoxEmitter::emit(MACGrid& grid, float dt){
    //1. 根据物理速率计算笨帧应该发射的粒子数
    float num_new_particles_float = m_particles_per_second * dt;
    int num_new_particles_int = static_cast<int>(num_new_particles_float);
    float fractional_part = num_new_particles_float - num_new_particles_int;
    //2. 用概率处理小数部分
    if (static_cast<float>(rand())/RAND_MAX < fractional_part){
        num_new_particles_int++;
    }
    //3. 循环创建新粒子
    for(int i = 0; i< num_new_particles_int; ++i){
        //长方体内生成随机位置
        float rand_x = static_cast<float>(rand()) / RAND_MAX;
        float rand_y = static_cast<float>(rand()) / RAND_MAX;
        float rand_z = static_cast<float>(rand()) / RAND_MAX;
        Vector3D particle_pos = m_min_corner + Vector3D(rand_x * (m_max_corner.x - m_min_corner.x),
                                                        rand_y * (m_max_corner.y - m_min_corner.y),
                                                        rand_z * (m_max_corner.z - m_min_corner.z));
        //创建一个新粒子
        Particles new_particle;
        new_particle.position = particle_pos;
        new_particle.velocity = m_initial_velocity;
        //将粒子添加到模拟系统中
        grid.particles().push_back(new_particle);
    }
}