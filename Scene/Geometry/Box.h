#pragma once
#include "SolidShape.h"
#include <algorithm>
#include <cmath>

class BoxShape:public SolidShape {
private:
    Point3f center;
    Point3f halfExtents;//长方体长宽高的一半

public:
    //构造函数: 传入中心点和半长宽高等几何信息
    BoxShape(const Point3f& c, const Point3f& h):center(c),halfExtents(h){}
    //SDF
    float signedDistance(const Point3f& position)const override{
        //计算点到中心的绝对距离, 并减去半长宽
        float dx=std::abs(position.x - center.x) - halfExtents.x;
        float dy=std::abs(position.y - center.y) - halfExtents.y;
        float dz=std::abs(position.z - center.z) - halfExtents.z;

        //外部距离(如果点在内部, 这一项为0)
        float out_dist=std::sqrt(std::max(dx, 0.0f) * std::max(dx, 0.0f) + 
                                std::max(dy, 0.0f) * std::max(dy, 0.0f) + 
                                std::max(dz, 0.0f) * std::max(dz, 0.0f));
        
        //内部距离(如果点在外部, 这一项为0)
        float in_dist = std::min(std::max({dx, dy, dz}), 0.0f);

        return out_dist + in_dist;
    }

    //实现法线计算(利用SDF的梯度, 这里使用中心差分法近似求导)
    Normal3f normal(const Point3f& position) const override {
        const float eps = 1e-4f; // 一个很小的偏移量用于求导
        
        // 分别对 x, y, z 方向求偏导数
        float nx = signedDistance(Point3f(position.x + eps, position.y, position.z)) - 
                   signedDistance(Point3f(position.x - eps, position.y, position.z));
        float ny = signedDistance(Point3f(position.x, position.y + eps, position.z)) - 
                   signedDistance(Point3f(position.x, position.y - eps, position.z));
        float nz = signedDistance(Point3f(position.x, position.y, position.z + eps)) - 
                   signedDistance(Point3f(position.x, position.y, position.z - eps));

        // 归一化法向量
        float length = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (length > 1e-6f) {
            return Normal3f(nx / length, ny / length, nz / length);
        }
        
        // 防退化处理（如果在完全正中心处等极端情况）
        return Normal3f(0.0f, 1.0f, 0.0f); 
    }
};