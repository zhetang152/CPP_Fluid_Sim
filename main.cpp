#include "repo.h"
void apply_pressure_gradient(MACGrid& grid, float dt){
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
void apply_pressure_gradient_FVM(MACGrid& grid, float dt) {
    const int nx = grid.getDimX();
    const int ny = grid.getDimY();
    const int nz = grid.getDimZ();
    const float scale = dt / (grid.getDx() * grid.density()(0,0,0)); // 假设密度恒定

    const auto& pressure = grid.pressure();
    const auto& u_area = grid.area_u();
    const auto& v_area = grid.area_v();
    const auto& w_area = grid.area_w();

    auto& u = grid.u();
    auto& v = grid.v();
    auto& w = grid.w();

    // 更新速度u
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 1; i < nx; ++i) {
                //面积分数 > 0, 就更新速度
                if (u_area(i, j, k) > 0) {
                    u(i, j, k) -= scale * (pressure(i, j, k) - pressure(i - 1, j, k));
                }
            }
        }
    }
    // 更新速度v
    for (int k = 0; k < nz; ++k) {
        for (int j = 1; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (v_area(i, j, k) > 0) {
                    v(i, j, k) -= scale * (pressure(i, j, k) - pressure(i, j - 1, k));
                }
            }
        }
    }
    // 更新速度w
    for (int k = 1; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (w_area(i, j, k) > 0) {
                    w(i, j, k) -= scale * (pressure(i, j, k) - pressure(i, j, k - 1));
                }
            }
        }
    }
}
int main(){
    //随机数种子
    std::srand(time(NULL));
    //1. 初始化
    const int resolution = 32;//分辨率
    MACGrid grid(resolution, resolution, resolution, 1.0f / static_cast<float>(resolution));
    std::unique_ptr<BoundaryCondition> boundary = std::make_unique<SolidBoundary>();
    std::vector<std::unique_ptr<ExternalForce>> forces;
    forces.emplace_back(std::make_unique<GravityForce>(Vector3D(0.0f, -9.81f, 0.0f))); //添加重力
    //模拟参数
    const float dt = 0.01f; //时间步长
    const float rho = 1.0f; //流体密度
    const int max_iterations = 200; //最大迭代次数
    const float tolerance = 1e-5f; //收敛容忍度
    #define USE_FVM_SOLVER 1
    //2. 设置场景
    SceneManager::createFishTank(grid, 0.8f); 
    #if USE_FVM_SOLVER
        Sphere solid_sphere(Vector3D(0.5, 0.5, 0.5), 0.15); 
        computeFractions(grid, solid_sphere, 4); 
    #endif
    //3. 主模拟循环
    for (int frame = 0; frame < 300; ++frame){
        std::cout << "Simulating Frame" << frame << std::endl;
        //advection
        Advector::advect_velocity(grid,dt);
        Advector::advect_particles(grid,solid_sphere,dt);
        //forces
        for (const auto& force : forces) {
            force->apply(grid, dt);
        }
        //A: 强制边界条件
        boundary->apply(grid);
        //B: 压力投影
        #if USE_FVM_SOLVER
            std::cout << "Using FVM Solver..." << std::endl;
            // B1: 计算FVM散度
            Grid<float> rhs = Solver::discrete_divergence_FVM(grid);
            // B2: 构建FVM矩阵A和预条件子
            Solver::SystemMatrix matrix(resolution, resolution, resolution);
            Solver::buildMatrixA_FVM(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid, dt);
            Solver::MIC0preconditioner_FVM(matrix.precon, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
            // B3: PCG_FVM求解压力
            Solver::PCG_FVM(grid.pressure(), rhs, matrix, grid, max_iterations, tolerance);
        #else
            //B1: 计算速度散度
            Grid<float> rhs = Solver::discrete_divergence(grid);
            //B2: 构建压力矩阵A和preconditioner
            Solver::SystemMatrix matrix(resolution, resolution, resolution);
            Solver::buildMatrixA(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes(), grid.getDx(), dt, rho);
            Solver::MIC0preconditioner(matrix.precon, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes());
            //B3: PCG求解压力
            Solver::PCG(grid.pressure(), rhs, matrix, grid.celltypes(), grid.getDx(), max_iterations, tolerance);
        #endif
        //C: 更新速度
        #if USE_FVM_SOLVER
            apply_pressure_gradient_FVM(grid, dt);
        #else
            apply_pressure_gradient(grid, dt);
        #endif
        //D: 再次强制边界条件
        boundary->apply(grid);
        //输出可视化文件
        std::stringstream ss;
        ss << "D:\\Computation\\FluidSim_result\\frame_" << std::setw(4) << std::setfill('0') << frame << ".obj";
        DataExporter::exportToObj(grid, ss.str());
    }
    std::cout << "Simulation completed." << std::endl;
    return 0;
}