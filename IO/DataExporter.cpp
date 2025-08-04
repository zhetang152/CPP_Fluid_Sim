#include "DataExporter.h"
#include<fstream>//用于文件操作
#include<iostream>

namespace DataExporter 
{
    bool exportToObj(const MACGrid& grid, const std::string& filepath){
        // 创建或打开文件流
        std::ofstream outfile(filepath);
        // 检查文件是否成功打开
        if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            return false;
        }
        const auto& particals = grid.particles();
        // 写入OBJ文件头
        outfile << "Partical data exported from FluidSim\n";
        outfile << "# Total particles: " << particals.size() <<"\n";
        //遍历粒子
        for (const auto& particle : particals) {
            outfile << "v " << particle.x << " " << particle.y << " " << particle.z << "\n";
        }
        //析构
        outfile.close();
        return true;
    }
} // namespace DataExporter
