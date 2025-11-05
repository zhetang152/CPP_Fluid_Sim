#include "Sphere.h"

Sphere::Sphere(const Point3f& center, float radius)
    : m_center(center), m_radius(radius) {}

float Sphere::signedDistance(const Point3f& position) const {
    //SDF核心: 查询点到球心的距离减去半径
    return (position - m_center).length() - m_radius;
}