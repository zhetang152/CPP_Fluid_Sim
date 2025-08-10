#include "SDFUtils.h"
#include "GridAndParticleSystem.h"
#include "D:\Computation\FluidSim\CPP_Sim\Scene\Geometry\SolidShape.h"
#include<iostream>

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