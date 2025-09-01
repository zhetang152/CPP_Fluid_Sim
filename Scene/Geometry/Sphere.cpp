#include "Sphere.hpp"

Sphere::Sphere(const Vector3D& center, float radius)
    : m_center(center), m_radius(radius) {}

float Sphere::signedDistance(const Vector3D& position) const {
    //SDF核心: 查询点到球心的距离减去半径
    return (position - m_center).length() - m_radius;
}