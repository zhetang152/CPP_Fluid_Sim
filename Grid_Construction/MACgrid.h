#pragma once
#include"vector3D.h"
#include"grid.h"

enum class CellType{
    FLUID,
    AIR,
    SOLID
};

class MACGrid{
private:
    int x,y,z;
    float dx;
    Grid<CellType> m_celltypes;
    Grid<float> m_u;
    Grid<float> m_v;
    Grid<float> m_w;
    Grid<float> m_pressure;
public:
    MACGrid(int tx, int ty, int tz, float initial_dx):
    x(tx),y(ty),z(tz),dx(initial_dx),
    m_celltypes(tx,ty,tz,CellType::AIR),
    m_pressure(tx,ty,tz),
    m_u(tx + 1, ty, tz),
    m_v(tx, ty + 1, tz),
    m_w(tx, ty, tz + 1)
    {
    }
    //提供访问
    Grid<float>& pressure() {return m_pressure;}
    //使用例: MAC.pressure()(i, j, k) = 10.f;
    const Grid<float>& pressure() const {return m_pressure;}
    Grid<float>& u() { return m_u; }
    const Grid<float>& u() const { return m_u; }
    Grid<float>& v() { return m_v; }
    const Grid<float>& v() const { return m_v; }
    Grid<float>& w() { return m_w; }
    const Grid<float>& w() const { return m_w; }
    Grid<CellType>& celltypes() { return m_celltypes; }
    const Grid<CellType>& celltypes() const { return m_celltypes; }
    float getDx()const{return dx;}
    Vector3D PositionOfPressure(int i, int j,int k) const{
        return Vector3D(
            (static_cast<float>(i)+0.5f)*dx,
            (static_cast<float>(j) + 0.5f) * dx,
            (static_cast<float>(k) + 0.5f) * dx
        );
    }
    Vector3D positionOfV(int i, int j, int k) const{
        return Vector3D(
            (static_cast<float>(i) + 0.5f) * dx,
            static_cast<float>(j) * dx,
            (static_cast<float>(k) + 0.5f) * dx
        );
    }
    Vector3D positionOfW(int i, int j, int k) const {
        return Vector3D(
            (static_cast<float>(i) + 0.5f) * dx,
            (static_cast<float>(j) + 0.5f) * dx,
            static_cast<float>(k) * dx
        );
    }
};