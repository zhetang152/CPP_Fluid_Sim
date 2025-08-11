#pragma once
#include "Grid_Construction\GridAndParticleSystem.h"

namespace SceneManager {
    /**
     * @brief 初始化为鱼缸模型
     * @param grid 要修改的MAC网格
     * @param fluidLevel 流体填充的高度 (0.0 to 1.0)
     */
    void createFishTank(MACGrid& grid, float fluidLevel = 0.5f);
}