#include"SolidShape.hpp"

class Sphere final : public SolidShape {
private:
    Vector3D m_center;
    float m_radius;
public:
    Sphere(const Vector3D& center, float radius);
    float signedDistance(const Vector3D& position) const override;
    // Normal
    Vector3D normal(const Vector3D& position) const override{
        Vector3D dir = position - m_center;
        return dir.normalize();
    };
};