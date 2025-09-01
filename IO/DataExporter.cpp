#include "DataExporter.hpp"
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
            outfile << "v " << particle.position.x << " " << particle.position.y << " " << particle.position.z << "\n";
        }
        //析构
        outfile.close();
        return true;
    }
    bool exportMeshToObj(const TriangleMesh& mesh, const std::string& filepath){
       std::ofstream outfile(filepath);
       if (!outfile.is_open()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            return false;
        }
        outfile << "# Mesh data exported from FluidSim\n";
        outfile << "# Vertices: " << mesh.vertices.size() << "\n";
        outfile << "# Faces: " << mesh.faces.size() / 3 << "\n";
        // 写入所有顶点
        for (const auto& v : mesh.vertices) {
            outfile << "v " << v.x << " " << v.y << " " << v.z << "\n";
        }
        // 写入所有法线
        for (const auto& n : mesh.normals) {
            outfile << "vn " << n.x << " " << n.y << " " << n.z << "\n";
        }
        // 写入所有面片
        for (size_t i = 0; i < mesh.faces.size(); i += 3) {
            outfile << "f " << mesh.faces[i] + 1 << " " 
                          << mesh.faces[i+1] + 1 << " " 
                          << mesh.faces[i+2] + 1 << "\n";
        }

        outfile.close();
        return true;
    }
}
