
#include"vecmath.h"
#include"grid.h"

enum class CellType{
    AIR,
    FLUID,
    SOLID
};
struct Particles{
    /* data */
    Point3f position;
    Vector3f velocity;
};
class MACGrid{
private:
    int x,y,z;
    Float dx; //网格大小
    Grid<CellType> m_celltypes; //version 1 : 使用枚举类型表示每个单元格的类型
    Grid<Float> m_u; //x方向速度分量
    Grid<Float> m_v; //y方向速度分量
    Grid<Float> m_w; //z方向速度分量
    Grid<Float> m_pressure; //压力场
    Grid<Float> m_density; //流体密度
    Grid<Vector3f> m_solidvelocity; //固体边界速度
    Grid<Float> m_liquid_phi; //液体的SDF
    Grid<Float> m_volume_fractions;//体积分数
    Grid<Float> m_area_u, m_area_v, m_area_w; //用于计算速度的开放面积分数
    std::vector<Particles> m_particles; //用于标记粒子位置的向量
    
    public:
    MACGrid(int tx, int ty, int tz, Float initial_dx):
    x(tx),y(ty),z(tz),dx(initial_dx),
    m_celltypes(tx,ty,tz,CellType::AIR),
    m_pressure(tx,ty,tz),
    m_solidvelocity(tx, ty, tz, Vector3f(0,0,0)),
    m_u(tx + 1, ty, tz),
    m_v(tx, ty + 1, tz),
    m_w(tx, ty, tz + 1),
    m_liquid_phi(tx, ty, tz, 0.0f),
    m_volume_fractions(tx, ty, tz, 0.0f),
    m_area_u(tx + 1, ty, tz, 0.0f),
    m_area_v(tx, ty + 1, tz, 0.0f),
    m_area_w(tx, ty, tz + 1, 0.0f),
    m_density(tx,ty,tz,0.0f)
    {}
    //提供访问
    int getDimX() const { return x; };
    int getDimY() const { return y; };
    int getDimZ() const { return z; };
    Grid<Float>& pressure() {return m_pressure;}
    //使用例: MAC.pressure()(i, j, k) = 10.f;
    const Grid<Float>& pressure() const {return m_pressure;}
    
    Grid<Vector3f>& solidvelocity() { return m_solidvelocity; }
    const Grid<Vector3f>& solidvelocity() const { return m_solidvelocity; }
    std::vector<Particles>& particles() { return m_particles; }
    const std::vector<Particles>& particles() const { return m_particles; }
    Grid<Float>& liquid_phi() {return m_liquid_phi;}
    const Grid<Float>& liquid_phi() const {return m_liquid_phi;}
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
    Float getDx()const{return dx;}
    Point3f PositionOfPressure(int i, int j,int k) const{
        return Point3f(
            (static_cast<Float>(i)+0.5f)*dx,
            (static_cast<Float>(j) + 0.5f) * dx,
            (static_cast<Float>(k) + 0.5f) * dx
        );
    }
    Point3f positionOfU(int i, int j, int k) const {
        return Point3f(
            static_cast<Float>(i) * dx,
            (static_cast<Float>(j) + 0.5f) * dx,
            (static_cast<Float>(k) + 0.5f) * dx
        );
    }
    Point3f positionOfV(int i, int j, int k) const{
        return Point3f(
            (static_cast<Float>(i) + 0.5f) * dx,
            static_cast<Float>(j) * dx,
            (static_cast<Float>(k) + 0.5f) * dx
        );
    }
    Point3f positionOfW(int i, int j, int k) const {
        return Point3f(
            (static_cast<Float>(i) + 0.5f) * dx,
            (static_cast<Float>(j) + 0.5f) * dx,
            static_cast<Float>(k) * dx
        );
    }
};