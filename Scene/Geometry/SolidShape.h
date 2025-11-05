

#include "Grid_Construction\vecmath.h"
class SolidShape {
public:
    virtual ~SolidShape() = default;
    //SDF
    virtual float signedDistance(const Point3f& position) const = 0;
    //Normal
    virtual Point3f normal(const Point3f& Point3fposition) const = 0;
};