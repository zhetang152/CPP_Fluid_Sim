#pragma once
#include "Force\ExternalForce.h"
#include "Grid_Construction\vector3D.h"
class GravityForce : public ExternalForce {
private:
    Vector3D m_gravity;
public:
    explicit GravityForce(const Vector3D& gravity);
    void apply(MACGrid& grid, float dt) override;
};