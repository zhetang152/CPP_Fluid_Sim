#include"SolidShape.h"
#include "Grid_Construction\math.h"

class Sphere final : public SolidShape {
private:
    Point3f m_center;
    Float m_radius;
public:
    Sphere(const Point3f& center, Float radius);
    float signedDistance(const Point3f& position) const override;
    // Normal
    Vector3f normal(const Point3f& position) const override{
        Vector3f dir = position - m_center;
        return dir.normalize();
    };
};