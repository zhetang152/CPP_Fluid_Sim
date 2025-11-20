#include "SceneManager.hpp"
#include<iostream>

namespace SceneManager {
    void createFishTank(MACGrid& grid, float fluidLevel){
        //基本参数
        const int nx = grid.getDimX();
        const int ny = grid.getDimY();
        const int nz = grid.getDimZ();
        Grid<CellType>& cellTypes = grid.celltypes();
        //获取固体速度
        Grid<Vector3f>& solidVelocity = grid.solidvelocity();
        //遍历所有网格
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    //检查单元格是否位于网格边界
                    if (i == 0 || i == nx - 1 || j == 0 || k == 0 || k == nz - 1) {
                        //设置为固体单元格
                        grid.celltypes()(i, j, k) = CellType::SOLID;
                        //设置固体速度为零
                        solidVelocity(i, j, k) = Vector3f(0.0f, 0.0f, 0.0f);
                    } else {
                        int fluidHeight = static_cast<int>(ny* fluidLevel);
                        if(j < fluidHeight){
                            cellTypes(i, j, k) = CellType::FLUID;
                            // 在每个被确定为 FLUID 的单元格中播种8个粒子
                            for (int p = 0; p < 8; ++p){
                                float jitter_x = (static_cast<float>(rand()) / RAND_MAX - 0.5f);
                                float jitter_y = (static_cast<float>(rand()) / RAND_MAX - 0.5f);
                                float jitter_z = (static_cast<float>(rand()) / RAND_MAX - 0.5f);
                                // 计算粒子位置，并添加到grid的粒子列表中
                                Point3f pos = grid.PositionOfPressure(i, j, k) + Vector3f(jitter_x, jitter_y, jitter_z) * grid.getDx();
                                Particles new_particle;
                                new_particle.position = pos;
                                new_particle.velocity = Vector3f(0.0f, 0.0f, 0.0f);
                                grid.particles().push_back(new_particle);
                            }
                        }else{
                            cellTypes(i, j, k) = CellType::AIR;
                        }
                    }
                }
            }
        }
        std::cout << "Scene 'FishTank' created with fluid level " << fluidLevel << "." << std::endl;
    };
}