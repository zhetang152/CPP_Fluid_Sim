#include "SDFUtils.hpp"
#include "GridAndParticleSystem.hpp"
#include "Scene\Geometry\SolidShape.hpp"
#include<iostream>
#include<vector>
#include<limits>

void computeFractions(MACGrid& grid, const SolidShape& solidShape, int supersample_level) {
    std::cout<< "Starting fraction computation with supersample level "<< supersample_level << std::endl;
    const float sub_dx = grid.getDx() / supersample_level;
    //计算体积分数
    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                int fluid_samples = 0;
                for (int sub_k = 0; sub_k < supersample_level; ++sub_k) {
                    for (int sub_j = 0; sub_j < supersample_level; ++sub_j) {
                        for (int sub_i = 0; sub_i < supersample_level; ++sub_i) {
                            float px = (i + (sub_i + 0.5f) / supersample_level) * grid.getDx();
                            float py = (j + (sub_j + 0.5f) / supersample_level) * grid.getDx();
                            float pz = (k + (sub_k + 0.5f) / supersample_level) * grid.getDx();
                            if (solidShape.signedDistance({px, py, pz}) > 0.0f) {
                                fluid_samples++;
                            }
                        }
                    }
                }
                float fraction = static_cast<float>(fluid_samples) / (supersample_level * supersample_level * supersample_level);
                grid.volumeFractions()(i, j, k) = fraction;
            }
        }
    }
    //计算面积分数
    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX() + 1; ++i) {
                int fluid_samples = 0;
                for (int sub_k = 0; sub_k < supersample_level; ++sub_k) {
                    for (int sub_j = 0; sub_j < supersample_level; ++sub_j) {
                        float px = static_cast<float>(i) * grid.getDx();
                        float py = (j + (sub_j + 0.5f) / supersample_level) * grid.getDx();
                        float pz = (k + (sub_k + 0.5f) / supersample_level) * grid.getDx();
                        if (solidShape.signedDistance({px, py, pz}) > 0.0f) {
                            fluid_samples ++;
                        }
                    }
                }
                float fraction = static_cast<float>(fluid_samples) / (supersample_level * supersample_level);
                grid.area_u()(i, j, k) = fraction;
            }
        }
    }
    // 计算v方向的面积分数
    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY() + 1; ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                int fluid_samples = 0;
                for (int sub_k = 0; sub_k < supersample_level; ++sub_k) {
                    for (int sub_i = 0; sub_i < supersample_level; ++sub_i) {
                        float px = (i + (sub_i + 0.5f) / supersample_level) * grid.getDx();
                        float py = static_cast<float>(j) * grid.getDx();
                        float pz = (k + (sub_k + 0.5f) / supersample_level) * grid.getDx();
                        if (solidShape.signedDistance({px, py, pz}) > 0.0f) {
                            fluid_samples++;
                        }
                    }
                }
                float fraction = static_cast<float>(fluid_samples) / (supersample_level * supersample_level);
                grid.area_v()(i, j, k) = fraction;
            }
        }
    }
    // 计算w方向的面积分数
    for (int k = 0; k < grid.getDimZ() + 1; ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                int fluid_samples = 0;
                for (int sub_j = 0; sub_j < supersample_level; ++sub_j) {
                    for (int sub_i = 0; sub_i < supersample_level; ++sub_i) {
                        float px = (i + (sub_i + 0.5f) / supersample_level) * grid.getDx();
                        float py = (j + (sub_j + 0.5f) / supersample_level) * grid.getDx();
                        float pz = static_cast<float>(k) * grid.getDx();
                        if (solidShape.signedDistance({px, py, pz}) > 0.0f) {
                            fluid_samples++;
                        }
                    }
                }
                float fraction = static_cast<float>(fluid_samples) / (supersample_level * supersample_level);
                grid.area_w()(i, j, k) = fraction;
            }
        }
    }
}

//粒子重构SDF的实现
void updateupdateLiquidSDFFromParticles(MACGrid& grid){
    auto& liquid_phi = grid.liquid_phi();
    const auto& particles = grid.particles();
    const float dx = grid.getDx();
    // 粒子的有效半径, 通常设置为比网格稍大, 以确保SDF连续
    const float particle_radius = dx * 1.5f; 
    // 1. 初始化SDF为一个很大的正数
    const float max_dist = 3 * dx;
    liquid_phi.fill(max_dist);
    //2. 遍历每个网格点
    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                Vector3D cell_center = grid.PositionOfPressure(i, j, k);
                float min_dist_sq = max_dist * max_dist;
                //3. 找到一定范围内的最近粒子
                for (const auto& p : particles) {
                    float dist_sq = (p.position - cell_center).length();
                    if (dist_sq < min_dist_sq) {
                        min_dist_sq = dist_sq;
                    }
                }
                //4. 计算SDF值: 最近点距离-粒子半径
                liquid_phi(i, j, k) = min_dist_sq - particle_radius;
            }
        }
    }
}
//根据SDF更新单元格类型的实现
void updateCellTypesFromSDF(MACGrid& grid) {
    auto& cell_types = grid.celltypes();
    const auto& liquid_phi = grid.liquid_phi();

    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                // 只更新非固体单元格
                if (cell_types(i, j, k) != CellType::SOLID) {
                    if (liquid_phi(i, j, k) <= 0.0f) {
                        cell_types(i, j, k) = CellType::FLUID;
                    } else {
                        cell_types(i, j, k) = CellType::AIR;
                    }
                }
            }
        }
    }
}