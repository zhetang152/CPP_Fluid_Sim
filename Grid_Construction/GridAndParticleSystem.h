#pragma once

#include <vector>

#include "grid.h"
#include "vecmath.h"

enum class CellType {
    AIR,
    FLUID,
    SOLID,
};

struct Particles {
    Point3f position;
    Vector3f velocity;
};

class MACGrid {
private:
    int x;
    int y;
    int z;
    Float dx;

    Grid<CellType> m_celltypes;
    Grid<Float> m_u;
    Grid<Float> m_v;
    Grid<Float> m_w;
    Grid<Float> m_pressure;
    Grid<Float> m_density;
    Grid<Vector3f> m_solidvelocity;
    Grid<Float> m_liquid_phi;
    Grid<Float> m_volume_fractions;
    Grid<Float> m_area_u;
    Grid<Float> m_area_v;
    Grid<Float> m_area_w;
    std::vector<Particles> m_particles;

public:
    MACGrid(int tx, int ty, int tz, Float initial_dx)
        : x(tx),
          y(ty),
          z(tz),
          dx(initial_dx),
          m_celltypes(tx, ty, tz, CellType::AIR),
          m_u(tx + 1, ty, tz),
          m_v(tx, ty + 1, tz),
          m_w(tx, ty, tz + 1),
          m_pressure(tx, ty, tz),
          m_density(tx, ty, tz, 0.0f),
          m_solidvelocity(tx, ty, tz, Vector3f(0, 0, 0)),
          m_liquid_phi(tx, ty, tz, 0.0f),
          m_volume_fractions(tx, ty, tz, 0.0f),
          m_area_u(tx + 1, ty, tz, 0.0f),
          m_area_v(tx, ty + 1, tz, 0.0f),
          m_area_w(tx, ty, tz + 1, 0.0f) {}

    int getDimX() const { return x; }
    int getDimY() const { return y; }
    int getDimZ() const { return z; }
    Float getDx() const { return dx; }

    Grid<Float>& pressure() { return m_pressure; }
    const Grid<Float>& pressure() const { return m_pressure; }

    Grid<Vector3f>& solidvelocity() { return m_solidvelocity; }
    const Grid<Vector3f>& solidvelocity() const { return m_solidvelocity; }

    std::vector<Particles>& particles() { return m_particles; }
    const std::vector<Particles>& particles() const { return m_particles; }

    Grid<Float>& liquid_phi() { return m_liquid_phi; }
    const Grid<Float>& liquid_phi() const { return m_liquid_phi; }

    Grid<Float>& volumeFractions() { return m_volume_fractions; }
    const Grid<Float>& volumeFractions() const { return m_volume_fractions; }

    Grid<Float>& area_u() { return m_area_u; }
    const Grid<Float>& area_u() const { return m_area_u; }

    Grid<Float>& area_v() { return m_area_v; }
    const Grid<Float>& area_v() const { return m_area_v; }

    Grid<Float>& area_w() { return m_area_w; }
    const Grid<Float>& area_w() const { return m_area_w; }

    Grid<Float>& u() { return m_u; }
    const Grid<Float>& u() const { return m_u; }

    Grid<Float>& v() { return m_v; }
    const Grid<Float>& v() const { return m_v; }

    Grid<Float>& w() { return m_w; }
    const Grid<Float>& w() const { return m_w; }

    Grid<Float>& density() { return m_density; }
    const Grid<Float>& density() const { return m_density; }

    Grid<CellType>& celltypes() { return m_celltypes; }
    const Grid<CellType>& celltypes() const { return m_celltypes; }

    Point3f PositionOfPressure(int i, int j, int k) const {
        return Point3f(
            (static_cast<Float>(i) + 0.5f) * dx,
            (static_cast<Float>(j) + 0.5f) * dx,
            (static_cast<Float>(k) + 0.5f) * dx);
    }

    Point3f positionOfU(int i, int j, int k) const {
        return Point3f(
            static_cast<Float>(i) * dx,
            (static_cast<Float>(j) + 0.5f) * dx,
            (static_cast<Float>(k) + 0.5f) * dx);
    }

    Point3f positionOfV(int i, int j, int k) const {
        return Point3f(
            (static_cast<Float>(i) + 0.5f) * dx,
            static_cast<Float>(j) * dx,
            (static_cast<Float>(k) + 0.5f) * dx);
    }

    Point3f positionOfW(int i, int j, int k) const {
        return Point3f(
            (static_cast<Float>(i) + 0.5f) * dx,
            (static_cast<Float>(j) + 0.5f) * dx,
            static_cast<Float>(k) * dx);
    }
};
