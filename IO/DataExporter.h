#pragma once
#include <string>
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\GridAndParticleSystem.h"
namespace DataExporter {
    /**
     * @brief 导出网格数据到.obj文件
     * @param grid 
     * @param filepath "D:\Computation\FluidSim_result\frame_0001.obj"
     * @return true / false
     */
    bool exportToObj(const MACGrid& grid, const std::string& filepath);
}