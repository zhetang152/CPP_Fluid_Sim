#pragma once

class MACGrid;
/**
 * 通用接口, 用于定义边界条件
 * 该接口可以被不同的边界条件实现类继承
 */
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;
    // 纯虚函数，任何继承此接口的派生类都必须自己实现 apply 函数
    virtual void apply(MACGrid& grid) = 0;
};