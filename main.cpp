#include<iostream>
#include<memory>
#include "D:\Computation\FluidSim\CPP_Sim\Grid_Construction\MACGrid.h"
#include "D:\Computation\FluidSim\CPP_Sim\Scene\SceneManager.h"
#include "D:\Computation\FluidSim\CPP_Sim\Boundary\SolidBoundary.h"
#include "D:\Computation\FluidSim\CPP_Sim\Solver\solver.h"

void apply_pressure_gardient(MACGrid& grid, float dt){
    const int nx = grid.getDimX();
    const int ny = grid.getDimY();
    const int nz = grid.getDimZ();
    const float scale = dt / grid.getDx();
    const auto& pressure = grid.pressure();
    const auto& celltypes = grid.celltypes();
    auto& u = grid.u();
    auto& v = grid.v();
    auto& w = grid.w();
    //更新速度u
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 1; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i-1, j, k) == CellType::FLUID) {
                    u(i, j, k) -= scale * (pressure(i, j, k) - pressure(i - 1, j, k));
                }
            }
        }
    }
    //更新速度v
    for (int k = 0; k < nz; ++k) {
        for (int j = 1; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i, j-1, k) == CellType::FLUID) {
                    v(i, j, k) -= scale * (pressure(i, j, k) - pressure(i, j - 1, k));
                }
            }
        }
    }
    //更新速度w
    for (int k = 1; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i, j, k-1) == CellType::FLUID) {
                    w(i, j, k) -= scale * (pressure(i, j, k) - pressure(i, j, k - 1));
                }
            }
        }
    }
}
int main(){
    //1. 初始化
    const int resolution = 64;//分辨率
    MACGrid grid(resolution, resolution, resolution, 1.0f / static_cast<float>(resolution));
    std::unique_ptr<BoundaryCondition> boundary = std::make_unique<SolidBoundary>();//?
    //模拟参数
    const float dt = 0.01f; //时间步长
    const float rho = 1.0f; //流体密度
    const int max_iterations = 200; //最大迭代次数
    const float tolerance = 1e-5f; //收敛容忍度
    //2. 设置场景
    SceneManager::createFishTank(grid, 0.5f); 
    grid.u()(resolution / 2, resolution / 2, resolution / 2) = 1.0f; //设置初始速度u
    grid.v()(resolution / 2, resolution / 2, resolution / 2) = -2.0f; //设置初始速度v
    //3. 主模拟循环
    for (int frame = 0; frame < 200; ++frame){
        std::cout << "Simulating Frame" << frame << std::endl;
        //A: 强制边界条件
        boundary->apply(grid);
        //B: 压力投影
        //B1: 计算速度散度
        Grid<float> rhs = Solver::discrete_divergence(grid);
        //B2: 构建压力矩阵A和preconditioner
        Solver::SystemMatrix matrix(resolution, resolution, resolution);
        Solver::buildMatrixA(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes(), grid.getDx(), dt, rho);
        Solver::MIC0preconditioner(matrix.precon, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes());
        //B3: PCG求解压力
        Solver::PCG(grid.pressure(), rhs, matrix, grid.celltypes(), grid.getDx(), max_iterations, tolerance);
        //C: 更新速度
        apply_pressure_gardient(grid, dt);
    }
    std::cout << "Simulation completed." << std::endl;
    return 0;
}