#include "SolidBoundary.h"

#include "Grid_Construction/GridAndParticleSystem.h"

void SolidBoundary::apply(MACGrid& grid) {
    const int nx = grid.getDimX();
    const int ny = grid.getDimY();
    const int nz = grid.getDimZ();

    Grid<CellType>& celltypes = grid.celltypes();
    Grid<Vector3f>& solidvelocity = grid.solidvelocity();
    Grid<float>& u = grid.u();
    Grid<float>& v = grid.v();
    Grid<float>& w = grid.w();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) != CellType::FLUID) {
                    continue;
                }

                if (i > 0 && celltypes(i - 1, j, k) == CellType::SOLID) {
                    u(i, j, k) = solidvelocity(i - 1, j, k).x;
                }
                if (i < nx - 1 && celltypes(i + 1, j, k) == CellType::SOLID) {
                    u(i + 1, j, k) = solidvelocity(i + 1, j, k).x;
                }

                if (j > 0 && celltypes(i, j - 1, k) == CellType::SOLID) {
                    v(i, j, k) = solidvelocity(i, j - 1, k).y;
                }
                if (j < ny - 1 && celltypes(i, j + 1, k) == CellType::SOLID) {
                    v(i, j + 1, k) = solidvelocity(i, j + 1, k).y;
                }

                if (k > 0 && celltypes(i, j, k - 1) == CellType::SOLID) {
                    w(i, j, k) = solidvelocity(i, j, k - 1).z;
                }
                if (k < nz - 1 && celltypes(i, j, k + 1) == CellType::SOLID) {
                    w(i, j, k + 1) = solidvelocity(i, j, k + 1).z;
                }
            }
        }
    }
}
