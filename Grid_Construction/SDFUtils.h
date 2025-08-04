#pragma once

class MACGrid;
class SolidShape;

/**
 * @brief 根据SDF计算并填充MACGrid中的体积分数和面积分数。
 * @param grid 要被写入分数的MACGrid对象。
 * @param solidShape 定义固体边界的SDF对象。
 * @param supersample_level 超采样的分辨率，例如传入 2 表示 2x2x2 的采样。
 */

 void computeFractions(MACGrid& grid, const SolidShape& solidShape, int supersample_level = 2);