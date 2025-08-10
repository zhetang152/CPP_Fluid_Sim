#pragma once
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\GridAndParticleSystem.h"

class ExternalForce {
public:
    virtual ~ExternalForce() = default;
    // apply函数接收grid和时间步长dt作为参数
    virtual void apply(MACGrid& grid, float dt) = 0;
};