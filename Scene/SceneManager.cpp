#include "SceneManager.h"
#include<iostream>

namespace SceneManager {
    void createFishTank(MACGrid& grid, float fluidLevel){
        //基本参数
        const int nx = grid.getDimX();
        const int ny = grid.getDimY();
        const int nz = grid.getDimZ();
        Grid<CellType>& cellTypes = grid.celltypes();
        //获取固体速度
        Grid<Vector3D>& solidVelocity = grid.solidvelocity();
        //遍历所有网格
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    //检查单元格是否位于网格边界
                    if (i == 0 || i == nx - 1 || j == 0 || k == 0 || k == nz - 1) {
                        //设置为固体单元格
                        grid.celltypes()(i, j, k) = CellType::SOLID;
                        //设置固体速度为零
                        solidVelocity(i, j, k) = Vector3D(0.0f, 0.0f, 0.0f);
                    } else {
                        int fluidHeight = static_cast<int>(ny* fluidLevel);
                        if(j < fluidHeight){
                            cellTypes(i, j, k) = CellType::FLUID;
                        }else{
                            cellTypes(i, j, k) = CellType::AIR;
                        }
                    }
                }
            }
            std::cout << "Scene 'FishTank' created with fluid level " << fluidLevel << "." << std::endl;
        }
    };
}