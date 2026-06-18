#include "SDFUtils.h"

#include <iostream>

#include "GridAndParticleSystem.h"
#include "Scene/Geometry/SolidShape.h"

void computeFractions(MACGrid& grid, const SolidShape& solidShape, int supersample_level) {
    std::cout << "Starting fraction computation with supersample level " << supersample_level << std::endl;
    const Float dx = grid.getDx();

    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                int fluid_samples = 0;
                for (int sub_k = 0; sub_k < supersample_level; ++sub_k) {
                    for (int sub_j = 0; sub_j < supersample_level; ++sub_j) {
                        for (int sub_i = 0; sub_i < supersample_level; ++sub_i) {
                            Float px = (i + (sub_i + 0.5f) / supersample_level) * dx;
                            Float py = (j + (sub_j + 0.5f) / supersample_level) * dx;
                            Float pz = (k + (sub_k + 0.5f) / supersample_level) * dx;
                            if (solidShape.signedDistance(Point3f(px, py, pz)) > 0.0f) {
                                ++fluid_samples;
                            }
                        }
                    }
                }
                grid.volumeFractions()(i, j, k) =
                    static_cast<Float>(fluid_samples) / (supersample_level * supersample_level * supersample_level);
            }
        }
    }

    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX() + 1; ++i) {
                int fluid_samples = 0;
                for (int sub_k = 0; sub_k < supersample_level; ++sub_k) {
                    for (int sub_j = 0; sub_j < supersample_level; ++sub_j) {
                        Float px = static_cast<Float>(i) * dx;
                        Float py = (j + (sub_j + 0.5f) / supersample_level) * dx;
                        Float pz = (k + (sub_k + 0.5f) / supersample_level) * dx;
                        if (solidShape.signedDistance(Point3f(px, py, pz)) > 0.0f) {
                            ++fluid_samples;
                        }
                    }
                }
                grid.area_u()(i, j, k) = static_cast<Float>(fluid_samples) / (supersample_level * supersample_level);
            }
        }
    }

    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY() + 1; ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                int fluid_samples = 0;
                for (int sub_k = 0; sub_k < supersample_level; ++sub_k) {
                    for (int sub_i = 0; sub_i < supersample_level; ++sub_i) {
                        Float px = (i + (sub_i + 0.5f) / supersample_level) * dx;
                        Float py = static_cast<Float>(j) * dx;
                        Float pz = (k + (sub_k + 0.5f) / supersample_level) * dx;
                        if (solidShape.signedDistance(Point3f(px, py, pz)) > 0.0f) {
                            ++fluid_samples;
                        }
                    }
                }
                grid.area_v()(i, j, k) = static_cast<Float>(fluid_samples) / (supersample_level * supersample_level);
            }
        }
    }

    for (int k = 0; k < grid.getDimZ() + 1; ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                int fluid_samples = 0;
                for (int sub_j = 0; sub_j < supersample_level; ++sub_j) {
                    for (int sub_i = 0; sub_i < supersample_level; ++sub_i) {
                        Float px = (i + (sub_i + 0.5f) / supersample_level) * dx;
                        Float py = (j + (sub_j + 0.5f) / supersample_level) * dx;
                        Float pz = static_cast<Float>(k) * dx;
                        if (solidShape.signedDistance(Point3f(px, py, pz)) > 0.0f) {
                            ++fluid_samples;
                        }
                    }
                }
                grid.area_w()(i, j, k) = static_cast<Float>(fluid_samples) / (supersample_level * supersample_level);
            }
        }
    }
}

void updateLiquidSDFFromParticles(MACGrid& grid) {
    auto& liquid_phi = grid.liquid_phi();
    const auto& particles = grid.particles();
    const Float dx = grid.getDx();
    const Float particle_radius = dx * 1.5f;
    const Float max_dist = 3 * dx;

    liquid_phi.fill(max_dist);

    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                Point3f cell_center = grid.PositionOfPressure(i, j, k);
                Float min_dist = max_dist;
                for (const auto& p : particles) {
                    Float dist = Length(p.position - cell_center);
                    if (dist < min_dist) {
                        min_dist = dist;
                    }
                }
                liquid_phi(i, j, k) = min_dist - particle_radius;
            }
        }
    }
}

void updateCellTypesFromSDF(MACGrid& grid) {
    auto& cell_types = grid.celltypes();
    const auto& liquid_phi = grid.liquid_phi();

    for (int k = 0; k < grid.getDimZ(); ++k) {
        for (int j = 0; j < grid.getDimY(); ++j) {
            for (int i = 0; i < grid.getDimX(); ++i) {
                if (cell_types(i, j, k) == CellType::SOLID) {
                    continue;
                }
                cell_types(i, j, k) = (liquid_phi(i, j, k) <= 0.0f) ? CellType::FLUID : CellType::AIR;
            }
        }
    }
}
