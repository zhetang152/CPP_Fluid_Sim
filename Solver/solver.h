#pragma once
#include"D:\Computation\FluidSim\CPP_Sim\Grid_Construction\grid.h"
#include"D:\Computation\FluidSim\CPP_Sim\Grid_Construction\MACgrid.h"
#include <iostream>
#include <cmath>
namespace Solver {
    //计算MACGrid的离散散度
    Grid<float> discrete_divergence(const MACGrid& grid);
    //基本的向量操作
    //dot product
    float dotProduct(const Grid<float>& a, const Grid<float>& b, const Grid<CellType>&cellTypes);
    //scalar product & vector addition
    void saxpy(Grid<float>& y, float a, const Grid<float>& x, const Grid<CellType>& cellTypes);
    //PCG
    //matrix-vector product
    void applyA(
    Grid<float>& result, 
    const Grid<float>& p,
    const Grid<CellType>& cellTypes,
    const Grid<float>& Adiag,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k
    );
    //buildmatrixA
    void buildMatrixA(
        Grid<float>& Adiag, 
        Grid<float>& Aplus_i, 
        Grid<float>& Aplus_j, 
        Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes,
        float dx, 
        float dt, 
        float rho
    );
    //preconditioner by MIC(0)
    void MIC0preconditioner(
        Grid<float>& precon,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes
    );
    //apply preconditioner
    void applyPreconditioner(
    Grid<float>& z,
    const Grid<float>& r,
    const Grid<float>& precon,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k,
    const Grid<CellType>& cellTypes
    );
    //PCG solver
    struct SystemMatrix {
        Grid<float> Adiag;
        Grid<float> Aplus_i;
        Grid<float> Aplus_j;
        Grid<float> Aplus_k;
        Grid<float> precon;
        SystemMatrix(int nx, int ny, int nz) :
            Adiag(nx, ny, nz, 0.0f),
            Aplus_i(nx, ny, nz, 0.0f),
            Aplus_j(nx, ny, nz, 0.0f),
            Aplus_k(nx, ny, nz, 0.0f),
            precon(nx, ny, nz, 0.0f)
        {}
    };
    void PCG(
        Grid<float>& p,
        const Grid<float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        float dx,
        int maxIterations,
        float tolerance
    );
}