#include "SceneManager.h"

#include <cstdlib>
#include <iostream>

namespace SceneManager {
    void createFishTank(MACGrid& grid, float fluidLevel) {
        const int nx = grid.getDimX();
        const int ny = grid.getDimY();
        const int nz = grid.getDimZ();
        const int fluidHeight = static_cast<int>(ny * fluidLevel);

        Grid<CellType>& cellTypes = grid.celltypes();
        Grid<Vector3f>& solidVelocity = grid.solidvelocity();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (i == 0 || i == nx - 1 || j == 0 || k == 0 || k == nz - 1) {
                        cellTypes(i, j, k) = CellType::SOLID;
                        solidVelocity(i, j, k) = Vector3f(0.0f, 0.0f, 0.0f);
                        continue;
                    }

                    if (j < fluidHeight) {
                        cellTypes(i, j, k) = CellType::FLUID;
                        for (int p = 0; p < 8; ++p) {
                            float jitter_x = static_cast<float>(rand()) / RAND_MAX - 0.5f;
                            float jitter_y = static_cast<float>(rand()) / RAND_MAX - 0.5f;
                            float jitter_z = static_cast<float>(rand()) / RAND_MAX - 0.5f;

                            Particles new_particle;
                            new_particle.position =
                                grid.PositionOfPressure(i, j, k) + Vector3f(jitter_x, jitter_y, jitter_z) * grid.getDx();
                            new_particle.velocity = Vector3f(0.0f, 0.0f, 0.0f);
                            grid.particles().push_back(new_particle);
                        }
                    } else {
                        cellTypes(i, j, k) = CellType::AIR;
                    }
                }
            }
        }

        std::cout << "Scene 'FishTank' created with fluid level " << fluidLevel << "." << std::endl;
    }
}
