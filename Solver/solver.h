#ifndef SOLVER_H
#define SOLVER_H
#include "Grid_Construction/grid.h"
#include "Grid_Construction/GridAndParticleSystem.h"
#include "Grid_Construction/math.h"
#include "Grid_Construction/vecmath.h"
#include <iostream>
#include <cmath>

namespace Solver {
    //基本的向量操作
    //dot product
    Float dotProduct(const Grid<Float>& a, const Grid<Float>& b, const Grid<CellType>&cellTypes);
    //scalar product & vector addition
    void saxpy(Grid<Float>& y, Float a, const Grid<Float>& x, const Grid<CellType>& cellTypes);
    //PCG
    //matrix-vector product
    void applyA(
    Grid<Float>& result, 
    const Grid<Float>& p,
    const Grid<CellType>& cellTypes,
    const Grid<Float>& Adiag,
    const Grid<Float>& Aplus_i,
    const Grid<Float>& Aplus_j,
    const Grid<Float>& Aplus_k
    );
    //计算MACGrid的离散散度
    Grid<Float> discrete_divergence(const MACGrid& grid);
    //buildmatrixA
    void buildMatrixA(
        Grid<Float>& Adiag, 
        Grid<Float>& Aplus_i, 
        Grid<Float>& Aplus_j, 
        Grid<Float>& Aplus_k,
        const Grid<CellType>& cellTypes,
        Float dx, 
        Float dt, 
        Float rho
    );
    Float dotProduct_FVM(const Grid<Float>& a, const Grid<Float>& b, const MACGrid& grid);
    void saxpy_FVM(Grid<Float>& y, Float a, const Grid<Float>& x, const MACGrid& grid);
    //[FVM & Variation] discrete divergence
    Grid<Float> discrete_divergence_FVM(const MACGrid& grid);
    //[FVM & Variation] build matrix A
    void buildMatrixA_FVM(
        Grid<Float>& Adiag, 
        Grid<Float>& Aplus_i, 
        Grid<Float>& Aplus_j, 
        Grid<Float>& Aplus_k,
        const MACGrid& grid,
        Float dt
    );
    void applyA_FVM(
        Grid<Float>& result, 
        const Grid<Float>& p,
        const MACGrid& grid,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k
    );
    //preconditioner by MIC(0)
    void MIC0preconditioner(
        Grid<Float>& precon,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const Grid<CellType>& cellTypes
    );
    //apply preconditioner
    void applyPreconditioner(
    Grid<Float>& z,
    const Grid<Float>& r,
    const Grid<Float>& precon,
    const Grid<Float>& Aplus_i,
    const Grid<Float>& Aplus_j,
    const Grid<Float>& Aplus_k,
    const Grid<CellType>& cellTypes
    );
    void MIC0preconditioner_FVM(
    Grid<Float>& precon,
    const Grid<Float>& Adiag,
    const Grid<Float>& Aplus_i,
    const Grid<Float>& Aplus_j,
    const Grid<Float>& Aplus_k,
    const MACGrid& grid // 传入 MACGrid 以获取 volumeFractions
);
    void applyPreconditioner_FVM(
        Grid<Float>& z,
        const Grid<Float>& r,
        const Grid<Float>& precon,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const MACGrid& grid
    );
    //PCG solver
    struct SystemMatrix {
        Grid<Float> Adiag;
        Grid<Float> Aplus_i;
        Grid<Float> Aplus_j;
        Grid<Float> Aplus_k;
        Grid<Float> precon;
        SystemMatrix(int nx, int ny, int nz) :
            Adiag(nx, ny, nz, 0.0f),
            Aplus_i(nx, ny, nz, 0.0f),
            Aplus_j(nx, ny, nz, 0.0f),
            Aplus_k(nx, ny, nz, 0.0f),
            precon(nx, ny, nz, 0.0f)
        {}
    };
    void PCG(
        Grid<Float>& p,
        const Grid<Float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        Float dx,
        int maxIterations,
        Float tolerance
    );
    void PCG_FVM(
        Grid<Float>& p,
        const Grid<Float>& b,
        const SystemMatrix& matrix,
        const MACGrid& grid,
        int maxIterations,
        Float tolerance
    );
}

#endif // SOLVER_H