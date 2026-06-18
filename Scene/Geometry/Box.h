#pragma once

#include <algorithm>
#include <cmath>

#include "SolidShape.h"

class BoxShape : public SolidShape {
private:
    Point3f center;
    Point3f halfExtents;

public:
    BoxShape(const Point3f& c, const Point3f& h) : center(c), halfExtents(h) {}

    float signedDistance(const Point3f& position) const override {
        float dx = std::abs(position.x - center.x) - halfExtents.x;
        float dy = std::abs(position.y - center.y) - halfExtents.y;
        float dz = std::abs(position.z - center.z) - halfExtents.z;

        float out_dist = std::sqrt(
            std::max(dx, 0.0f) * std::max(dx, 0.0f) +
            std::max(dy, 0.0f) * std::max(dy, 0.0f) +
            std::max(dz, 0.0f) * std::max(dz, 0.0f));
        float in_dist = std::min(std::max({dx, dy, dz}), 0.0f);
        return out_dist + in_dist;
    }

    Normal3f normal(const Point3f& position) const override {
        const float eps = 1e-4f;
        float nx = signedDistance(Point3f(position.x + eps, position.y, position.z)) -
                   signedDistance(Point3f(position.x - eps, position.y, position.z));
        float ny = signedDistance(Point3f(position.x, position.y + eps, position.z)) -
                   signedDistance(Point3f(position.x, position.y - eps, position.z));
        float nz = signedDistance(Point3f(position.x, position.y, position.z + eps)) -
                   signedDistance(Point3f(position.x, position.y, position.z - eps));

        float length = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (length > 1e-6f) {
            return Normal3f(nx / length, ny / length, nz / length);
        }
        return Normal3f(0.0f, 1.0f, 0.0f);
    }
};
