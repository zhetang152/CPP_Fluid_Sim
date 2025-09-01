#include"SolidBoundary.hpp"
#include "Grid_Construction\GridAndParticleSystem.hpp"

void SolidBoundary::apply(MACGrid& grid){
    const int nx = grid.getDimX();
    const int ny = grid.getDimY();
    const int nz = grid.getDimZ();
    //获取内部数据网格的引用
    Grid<CellType>& celltypes = grid.celltypes();
    Grid<Vector3D>& solidvelocity = grid.solidvelocity();
    Grid<float>& u = grid.u();
    Grid<float>& v = grid.v();
    Grid<float>& w = grid.w();
    //遍历所有单元格, 如果单元格是FLUID, 则检查neighborhood
    for(int k = 0;k<nz;++k){
        for (int j = 0;j<ny;++j){
            for (int i = 0;i<nx;++i){
                if (celltypes(i,j,k)==CellType::FLUID){
                    //x方向速度
                    if(i > 0 && celltypes(i - 1, j, k) == CellType::SOLID){
                        u(i, j, k) = solidvelocity(i-1,j,k).x;
                    }
                    if (i<nx-1 && celltypes(i + 1, j, k) == CellType::SOLID){
                        u(i + 1, j, k) = solidvelocity(i + 1, j, k).x;
                    }
                    //y方向速度
                    if (j > 0 && celltypes(i, j - 1, k) == CellType::SOLID) {
                        v(i, j, k) = solidvelocity(i, j - 1, k).y;
                    }
                    if (j < ny - 1 && celltypes(i, j + 1, k) == CellType::SOLID) {
                        v(i, j + 1, k) = solidvelocity(i, j + 1, k).y;
                    }
                    //z方向速度
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
}