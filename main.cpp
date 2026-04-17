
#include "index.h"

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
                    Float p_center = pressure(i, j, k);
                    Float p_left = pressure(i - 1, j, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_left;
                    if (celltypes(i - 1, j, k) == CellType::SOLID) p_left = p_center;
                    u(i, j, k) -= scale * (p_center - p_left);
                }
            }
        }
    }
    //更新速度v
    for (int k = 0; k < nz; ++k) {
        for (int j = 1; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i, j-1, k) == CellType::FLUID) {
                    Float p_center = pressure(i, j, k);
                    Float p_bottom = pressure(i, j - 1, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_bottom;
                    if (celltypes(i, j - 1, k) == CellType::SOLID) p_bottom = p_center;
                    v(i, j, k) -= scale * (p_center - p_bottom);
                }
            }
        }
    }
    //更新速度w
    for (int k = 1; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i, j, k-1) == CellType::FLUID) {
                   Float p_center = pressure(i, j, k);
                    Float p_back = pressure(i, j, k - 1);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_back;
                    if (celltypes(i, j, k - 1) == CellType::SOLID) p_back = p_center;
                    w(i, j, k) -= scale * (p_center - p_back);
                }
            }
        }
    }
}
void apply_pressure_gradient_FVM(MACGrid& grid, float dt, float rho) {
    const int nx = grid.getDimX();
    const int ny = grid.getDimY();
    const int nz = grid.getDimZ();
    const Float dx = grid.getDx();
    const auto& pressure = grid.pressure();
    const auto& u_area = grid.area_u();
    const auto& v_area = grid.area_v();
    const auto& w_area = grid.area_w();
    const auto& density = grid.density();
    const auto& celltypes = grid.celltypes();
    auto& u = grid.u();
    auto& v = grid.v();
    auto& w = grid.w();

    // 更新速度u
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 1; i < nx; ++i) {
                //面积分数 > 0, 就更新速度
                if (u_area(i, j, k) > 0) {
                    Float rho_face = 0.5f * (density(i, j, k) + density(i - 1, j, k));
                    if (rho_face < 1e-6f) rho_face = rho;
                    Float scale = dt / (dx * rho_face);\
                    Float p_center = pressure(i, j, k);
                    Float p_left = pressure(i - 1, j, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_left;
                    if (celltypes(i - 1, j, k) == CellType::SOLID) p_left = p_center;
                    u(i, j, k) -= scale * (p_center - p_left);
                }
            }
        }
    }
    // 更新速度v
    for (int k = 0; k < nz; ++k) {
        for (int j = 1; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (v_area(i, j, k) > 0) {
                    Float rho_face = 0.5f * (density(i, j, k) + density(i, j - 1, k));
                    if (rho_face < 1e-6f) rho_face = rho;
                    Float scale = dt / (dx * rho_face);
                    Float p_center = pressure(i, j, k);
                    Float p_bottom = pressure(i, j - 1, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_bottom;
                    if (celltypes(i, j - 1, k) == CellType::SOLID) p_bottom = p_center;
                    v(i, j, k) -= scale * (p_center - p_bottom);
                }
            }
        }
    }
    // 更新速度w
    for (int k = 1; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (w_area(i, j, k) > 0) {
                    Float rho_face = 0.5f * (density(i, j, k) + density(i, j, k - 1));
                    if (rho_face < 1e-6f) rho_face = rho;
                    Float scale = dt / (dx * rho_face);
                    Float p_center = pressure(i, j, k);
                    Float p_back = pressure(i, j, k - 1);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_back;
                    if (celltypes(i, j, k - 1) == CellType::SOLID) p_back = p_center;
                    w(i, j, k) -= scale * (p_center - p_back);
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
    const float dx = 1.0f / static_cast<float>(resolution);
    MACGrid grid(resolution, resolution, resolution, dx);
    std::unique_ptr<BoundaryCondition> boundary = std::make_unique<SolidBoundary>();
    std::vector<std::unique_ptr<ExternalForce>> forces;
    forces.emplace_back(std::make_unique<GravityForce>(Vector3f(0.0f, -9.81f, 0.0f))); //添加重力
    
    //创建粒子发射器
    std::vector<std::unique_ptr<ParticleEmitter>> emitters;
    emitters.emplace_back(std::make_unique<BoxEmitter>(
        Point3f(0.4,0.8,0.4),//发射区域的最小角
        Point3f(0.6,1.0,0.6),//发射区域的最大角
        2000.0f,
        Vector3f(0.0,-1.0,0.0)//初始速度向下
    ));
    //模拟参数
    const float dt = 0.01f; //时间步长
    const float rho = 1.0f; //流体密度
    const int max_iterations = 200; //最大迭代次数
    const float tolerance = 1e-5f; //收敛容忍度
    const float flip_alpha = 0.05f; // PIC/FLIP混合系数

    #define USE_FVM_SOLVER 1

    //创建并初始化marching cubes实例
    MarchingCubes mc(resolution, resolution, resolution);
    mc.set_int_data();
    mc.init_all();
    //2. 设置场景
    SceneManager::createFishTank(grid, 0.2f);
    BoxShape solid_box(Point3f(0.5f, 0.5f, 0.5f), Point3f(0.3f, 0.2f, 0.15f));

    #if USE_FVM_SOLVER
        computeFractions(grid, solid_box, 4); 
    #endif

    //3. 主模拟循环
    for(int frame = 0; frame < 300; ++frame){
        std::cout << "Simulating Frame" << frame << std::endl;
        //FLIP
        //发射新粒子
        for(const auto& emitter : emitters){
            emitter->emit(grid, dt);
        }
        //P2G
        FLIPSolver::ParticleToGrid(grid);
        //创建副本,计算变化量
        Grid<float> u_old = grid.u();
        Grid<float> v_old = grid.v();
        Grid<float> w_old = grid.w();
        //在网格上进行计算
        //  加外力
        for(const auto& force : forces){
            force->apply(grid, dt);
        }
        //  强制边界条件
        boundary->apply(grid);
        //  压力投影求解
        #if USE_FVM_SOLVER
            std::cout << "Using FVM Solver..." << std::endl;
            // B1: 计算FVM散度
            Grid<float> rhs = Solver::discrete_divergence_FVM(grid);
            // B2: 构建FVM矩阵A和预条件子
            Solver::SystemMatrix matrix(resolution, resolution, resolution);
            Solver::buildMatrixA_FVM(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid, grid.celltypes(), dt, rho);
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
            apply_pressure_gradient_FVM(grid, dt,rho);
        #else
            apply_pressure_gradient(grid, dt);
        #endif
        // ==========================================
        // STEP 2 新增：强制速度截断防爆逻辑
        // 在更新完速度后，检查并限制异常极大的速度
        // ==========================================
        float max_allowed_vel = 20.0f;
        auto& u_vel = grid.u();
        auto& v_vel = grid.v();
        auto& w_vel = grid.w();

        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    if (std::abs(u_vel(i, j, k)) > max_allowed_vel) {
                        u_vel(i, j, k) = (u_vel(i, j, k) > 0) ? max_allowed_vel : -max_allowed_vel;
                    }
                }
            }
        }
        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    if (std::abs(v_vel(i, j, k)) > max_allowed_vel) {
                        v_vel(i, j, k) = (v_vel(i, j, k) > 0) ? max_allowed_vel : -max_allowed_vel;
                    }
                }
            }
        }
        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    if (std::abs(w_vel(i, j, k)) > max_allowed_vel) {
                        w_vel(i, j, k) = (w_vel(i, j, k) > 0) ? max_allowed_vel : -max_allowed_vel;
                    }
                }
            }
        }
        // ==========================================
        //D: 再次强制边界条件
        boundary->apply(grid);
        //G2P
        FLIPSolver::GridToParticle(grid, u_old, v_old, w_old, flip_alpha);
        //平流粒子
        Advector::advect_particles(grid,solid_box, dt);
        //水平集更新
        updateLiquidSDFFromParticles(grid);
        updateCellTypesFromSDF(grid);
        //粒子管理
        const int min_particles_per_cell = 3;
        const int max_particles_per_cell = 12;
        // 创建一个网格来统计每个单元格的粒子数
        Grid<int> particle_counts(grid.getDimX(), grid.getDimY(), grid.getDimZ(), 0);
        for (const auto& p : grid.particles()) {
            int i = static_cast<int>(p.position.x / grid.getDx());
            int j = static_cast<int>(p.position.y / grid.getDx());
            int k = static_cast<int>(p.position.z / grid.getDx());
            // 安全边界检查
            if (i >= 0 && i < grid.getDimX() && j >= 0 && j < grid.getDimY() && k >= 0 && k < grid.getDimZ()) {
                particle_counts(i, j, k)++;
            }
        }
        // 2. 提前计算好每个单元格"必须"删除的粒子数量
        Grid<int> particles_to_drop(grid.getDimX(), grid.getDimY(), grid.getDimZ(), 0);
        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    if (grid.celltypes()(i, j, k) == CellType::FLUID && particle_counts(i, j, k) > max_particles_per_cell) {
                        particles_to_drop(i, j, k) = particle_counts(i, j, k) - max_particles_per_cell;
                    }
                }
            }
        }
        //删除拥堵区域的多余粒子
        auto& particles = grid.particles();
        particles.erase(std::remove_if(particles.begin(), particles.end(), 
            [&](const Particles& p) {
                int i = static_cast<int>(p.position.x / grid.getDx());
                int j = static_cast<int>(p.position.y / grid.getDx());
                int k = static_cast<int>(p.position.z / grid.getDx());

                if (i < 0 || i >= grid.getDimX() || j < 0 || j >= grid.getDimY() || k < 0 || k >= grid.getDimZ()) {
                    return true; // 删除跑到边界外的粒子
                }

                // 只在流体单元格中考虑删除

                if (grid.celltypes()(i, j, k) == CellType::FLUID) {
                    int& c = particle_counts(i, j, k);    // 该单元格目前还剩下几个粒子等待评估
                    int& d = particles_to_drop(i, j, k);  // 该单元格目前还需要删除几个粒子

                    if (d > 0 && c > 0) { // 如果还需要删，并且还有候选的
                        // 计算当前这颗粒子的删除概率
                        float prob = static_cast<float>(d) / static_cast<float>(c);
                        c--; // 不管删不删，这个粒子评估过了，剩余候选数-1
                        
                        if (static_cast<float>(rand()) / RAND_MAX < prob) {
                            d--; // 命中概率，所需删除数量-1
                            return true; // 返回true告诉vector删除此粒子
                        }
                    } else if (c > 0) {
                        c--; // 目标数量已经达标，后续粒子直接放行，但仍要减少剩余候选数
                    }
                }
                return false;
            }), particles.end());
        /*
        //在稀疏区域补充新粒子 (重播种)
        std::vector<Particles> new_particles;
        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    // 只在流体单元格中考虑补充
                    if (grid.celltypes()(i, j, k) == CellType::FLUID && 
                        grid.liquid_phi()(i, j, k) < -0.5f * grid.getDx() && 
                        particle_counts(i, j, k) < min_particles_per_cell) {
                        int num_to_add = min_particles_per_cell - particle_counts(i, j, k);
                        for (int n = 0; n < num_to_add; ++n) {
                            // 使用抖动随机创建新粒子
                            float jitter_x = (static_cast<float>(rand()) / RAND_MAX);
                            float jitter_y = (static_cast<float>(rand()) / RAND_MAX);
                            float jitter_z = (static_cast<float>(rand()) / RAND_MAX);
                            
                            Point3f pos = Point3f((i + jitter_x) * grid.getDx(), (j + jitter_y) * grid.getDx(), (k + jitter_z) * grid.getDx());

                            // 创建新粒子并从网格插值速度
                            Particles new_p;
                            new_p.position = pos;
                            new_p.velocity = Advector::get_velocity_at(grid, pos); // 从当前网格速度插值
                            new_particles.push_back(new_p);
                        }
                    }
                }
            }
        }
        // 将所有新创建的粒子一次性加入主列表
        particles.insert(particles.end(), new_particles.begin(), new_particles.end());
        */
        //A. 从SDF生成用于渲染的三角网格
        mc.reset();
        //B. 将SDF数据传输到MarchingCubes
        for(int k = 0; k < grid.getDimZ(); ++k){
            for(int j = 0; j < grid.getDimY(); ++j){
                for(int i = 0; i < grid.getDimX(); ++i){
                    mc.set_data(grid.liquid_phi()(i,j,k), i, j, k);
                }
            }
        }
        //C. 运行算法, 生成网格
        mc.set_method(false);
        mc.run(0.0f);
        //D. 将结果进行转换
        TriangleMesh render_mesh;
        render_mesh.vertices.reserve(mc.nverts());
        render_mesh.normals.reserve(mc.nverts());
        render_mesh.faces.reserve(mc.ntrigs() * 3);

        for (int i = 0; i < mc.nverts(); ++i) {
            Vertex* v = mc.vert(i);
            render_mesh.vertices.emplace_back(v->x * dx, v->y * dx, v->z * dx);
            render_mesh.normals.emplace_back(v->nx, v->ny, v->nz);
        }
        for (int i = 0; i < mc.ntrigs(); ++i) {
            Triangle* t = mc.trig(i);
            render_mesh.faces.push_back(t->v1);
            render_mesh.faces.push_back(t->v2);
            render_mesh.faces.push_back(t->v3);
        }
        //end
        //输出可视化文件
        std::stringstream ss;
        ss << "D:\\Code\\XLAB\\Fluid\\result\\frame_" << std::setw(4) << std::setfill('0') << frame << ".obj";
        //若要看粒子: DataExporter::exportToObj(grid, ss.str());
        //查看表面
        DataExporter::exportMeshToObj(render_mesh, ss.str());
    }
    std::cout << "Simulation completed." << std::endl;
    return 0;
}