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

/**
 * @brief 根据粒子位置重新计算液体SDF
 * @param grid 包含粒子和待更新m_liquid_phi的MACGrid
 * @param particle_radius 粒子的有效半径, 用于计算距离 
 */
void updateLiquidSDFFromParticles(MACGrid& grid);

/**
 * @brief 根据液体SDF更新单元格类型(FLUID/AIR)
 * @param grid 包含m_liquid_phi和待更新m_celltypes的MACGrid
 */
void updateCellTypesFromSDF(MACGrid& grid);