#pragma once
#include "Grid_Construction\vector3D.hpp"
class SolidShape {
public:
    virtual ~SolidShape() = default;
    //SDF
    virtual float signedDistance(const Vector3D& position) const = 0;
    //Normal
    virtual Vector3D normal(const Vector3D& position) const = 0;
};