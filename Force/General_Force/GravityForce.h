#pragma once
#include "Force\ExternalForce.h"
#include "Grid_Construction\vecmath.h"
class GravityForce : public ExternalForce {
private:
    Vector3f m_gravity;
public:
    explicit GravityForce(const Vector3f& gravity);
    void apply(MACGrid& grid, float dt) override;
};