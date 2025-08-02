#pragma once
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\vector3D.h"
class SolidShape {
public:
    virtual ~SolidShape() = default;
    //SDF
    virtual float signedDistance(const Vector3D& position) const = 0;
};