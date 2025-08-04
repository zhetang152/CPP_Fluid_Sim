#pragma once
#include "D:\Computation\FluidSim\CPP_Sim\Force\ExternalForce.h"
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\vector3D.h"
class GravityForce : public ExternalForce {
private:
    Vector3D m_gravity;
public:
    explicit GravityForce(const Vector3D& gravity);
    void apply(MACGrid& grid, float dt) override;
};