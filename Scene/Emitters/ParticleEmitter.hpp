#pragma once
#include"Grid_Construction\GridAndParticleSystem.hpp"
class ParticleEmitter {
public: 
    virtual ~ParticleEmitter()=default;
    /**
     * @brief 每个时间步都被调用, 用于发射粒子
     * @param grid 要将粒子发射到其中的MACGrid
     * @param dt 时间步长
    */
    virtual void emit(MACGrid& grid, float dt) = 0;
};