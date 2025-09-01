#pragma once
#include"BoundaryCondition.hpp"
/**
 * @class SolidBoundary
 * @brief 处理固体墙壁的速度边界条件
 */
class SolidBoundary : public BoundaryCondition {
public:
    SolidBoundary() = default;
    //override?
    void apply(MACGrid& grid) override;
};