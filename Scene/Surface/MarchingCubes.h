#pragma once
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\grid.h"
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\vector3D.h"
#include <vector>

// 用于存储最终三角网格的数据结构
struct TriangleMesh {
    std::vector<Vector3D> vertices;
    std::vector<Vector3D> normals; //用于平滑着色
    std::vector<int> faces;       //每3个整数代表一个三角形面
};
namespace MarchingCubes {
    /**
     * @brief 从SDF网格中提取等值面(isosurface)
     * @param sdf_grid 包含SDF值的网格
     * @param dx 网格单元尺寸
     * @param iso_level 要提取的等值面的值
     * @return 一个包含顶点和面片列表的TriangleMesh对象
     */
    TriangleMesh extractSurface(const Grid<float>& sdf_grid, float dx, float iso_level = 0.0f);
}