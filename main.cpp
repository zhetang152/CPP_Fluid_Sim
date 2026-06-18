#include "index.h"

#include <algorithm>
#include <cmath>
#include <filesystem>

void apply_pressure_gradient(MACGrid& grid, float dt) {
    const int nx = grid.getDimX();
    const int ny = grid.getDimY();
    const int nz = grid.getDimZ();
    const float scale = dt / grid.getDx();
    const auto& pressure = grid.pressure();
    const auto& celltypes = grid.celltypes();
    auto& u = grid.u();
    auto& v = grid.v();
    auto& w = grid.w();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 1; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i - 1, j, k) == CellType::FLUID) {
                    Float p_center = pressure(i, j, k);
                    Float p_left = pressure(i - 1, j, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_left;
                    if (celltypes(i - 1, j, k) == CellType::SOLID) p_left = p_center;
                    u(i, j, k) -= scale * (p_center - p_left);
                }
            }
        }
    }

    for (int k = 0; k < nz; ++k) {
        for (int j = 1; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i, j - 1, k) == CellType::FLUID) {
                    Float p_center = pressure(i, j, k);
                    Float p_bottom = pressure(i, j - 1, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_bottom;
                    if (celltypes(i, j - 1, k) == CellType::SOLID) p_bottom = p_center;
                    v(i, j, k) -= scale * (p_center - p_bottom);
                }
            }
        }
    }

    for (int k = 1; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                if (celltypes(i, j, k) == CellType::FLUID || celltypes(i, j, k - 1) == CellType::FLUID) {
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

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 1; i < nx; ++i) {
                if (u_area(i, j, k) > 0) {
                    Float rho_face = 0.5f * (density(i, j, k) + density(i - 1, j, k));
                    if (rho_face < 1e-6f) rho_face = rho;
                    Float scale = dt / (dx * rho_face);
                    Float p_center = pressure(i, j, k);
                    Float p_left = pressure(i - 1, j, k);
                    if (celltypes(i, j, k) == CellType::SOLID) p_center = p_left;
                    if (celltypes(i - 1, j, k) == CellType::SOLID) p_left = p_center;
                    u(i, j, k) -= scale * (p_center - p_left);
                }
            }
        }
    }

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

void clampGridValues(Grid<float>& values, float limit) {
    for (int k = 0; k < values.getDepth(); ++k) {
        for (int j = 0; j < values.getHeight(); ++j) {
            for (int i = 0; i < values.getWidth(); ++i) {
                values(i, j, k) = std::clamp(values(i, j, k), -limit, limit);
            }
        }
    }
}

int main() {
    std::srand(static_cast<unsigned int>(time(nullptr)));

    constexpr int resolution = 32;
    const float dx = 1.0f / static_cast<float>(resolution);
    constexpr float dt = 0.0025f;
    constexpr float rho = 1.0f;
    constexpr int max_iterations = 200;
    constexpr float tolerance = 1e-5f;
    constexpr float flip_alpha = 0.10f;
    constexpr bool use_fvm_solver = true;

    MACGrid grid(resolution, resolution, resolution, dx);
    std::unique_ptr<BoundaryCondition> boundary = std::make_unique<SolidBoundary>();

    std::vector<std::unique_ptr<ExternalForce>> forces;
    forces.emplace_back(std::make_unique<GravityForce>(Vector3f(0.0f, -9.81f, 0.0f)));

    std::vector<std::unique_ptr<ParticleEmitter>> emitters;
    emitters.emplace_back(std::make_unique<BoxEmitter>(
        Point3f(0.4f, 0.8f, 0.4f),
        Point3f(0.6f, 1.0f, 0.6f),
        2000.0f,
        Vector3f(0.0f, -1.0f, 0.0f)));

    MarchingCubes mc(resolution, resolution, resolution);
    mc.set_int_data();
    mc.init_all();

    SceneManager::createFishTank(grid, 0.2f);
    BoxShape solid_box(Point3f(0.5f, 0.5f, 0.5f), Point3f(0.3f, 0.2f, 0.15f));

    if constexpr (use_fvm_solver) {
        computeFractions(grid, solid_box, 4);
    }

    const std::filesystem::path output_dir("D:/Code/XLAB/Fluid/result");
    std::filesystem::create_directories(output_dir);

    for (int frame = 0; frame < 300; ++frame) {
        std::cout << "Simulating Frame " << frame << std::endl;

        for (const auto& emitter : emitters) {
            emitter->emit(grid, dt);
        }

        FLIPSolver::ParticleToGrid(grid);

        Grid<float> u_old = grid.u();
        Grid<float> v_old = grid.v();
        Grid<float> w_old = grid.w();

        for (const auto& force : forces) {
            force->apply(grid, dt);
        }

        boundary->apply(grid);
        grid.pressure().fill(0.0f);

        if constexpr (use_fvm_solver) {
            std::cout << "Using FVM Solver..." << std::endl;
            Grid<float> rhs = Solver::discrete_divergence_FVM(grid);
            Solver::SystemMatrix matrix(resolution, resolution, resolution);
            Solver::buildMatrixA_FVM(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid, grid.celltypes(), dt, rho);
            Solver::MIC0preconditioner_FVM(matrix.precon, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
            Solver::PCG_FVM(grid.pressure(), rhs, matrix, grid, max_iterations, tolerance);
            apply_pressure_gradient_FVM(grid, dt, rho);
        } else {
            Grid<float> rhs = Solver::discrete_divergence(grid);
            Solver::SystemMatrix matrix(resolution, resolution, resolution);
            Solver::buildMatrixA(matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes(), grid.getDx(), dt, rho);
            Solver::MIC0preconditioner(matrix.precon, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid.celltypes());
            Solver::PCG(grid.pressure(), rhs, matrix, grid.celltypes(), grid.getDx(), max_iterations, tolerance);
            apply_pressure_gradient(grid, dt);
        }

        constexpr float max_allowed_vel = 20.0f;
        clampGridValues(grid.u(), max_allowed_vel);
        clampGridValues(grid.v(), max_allowed_vel);
        clampGridValues(grid.w(), max_allowed_vel);

        boundary->apply(grid);

        FLIPSolver::GridToParticle(grid, u_old, v_old, w_old, flip_alpha);
        Advector::advect_particles(grid, solid_box, dt);

        const float margin = 0.001f;
        const float max_x = grid.getDimX() * grid.getDx() - margin;
        const float max_y = grid.getDimY() * grid.getDx() - margin;
        const float max_z = grid.getDimZ() * grid.getDx() - margin;
        const float max_speed = 15.0f;

        for (auto& p : grid.particles()) {
            p.position.x = std::clamp(p.position.x, margin, max_x);
            p.position.y = std::clamp(p.position.y, margin, max_y);
            p.position.z = std::clamp(p.position.z, margin, max_z);

            float speed_sq = p.velocity.x * p.velocity.x +
                             p.velocity.y * p.velocity.y +
                             p.velocity.z * p.velocity.z;
            if (speed_sq > max_speed * max_speed) {
                float speed = std::sqrt(speed_sq);
                p.velocity = p.velocity * (max_speed / speed);
            }
        }

        updateLiquidSDFFromParticles(grid);
        updateCellTypesFromSDF(grid);

        const int max_particles_per_cell = 12;
        Grid<int> particle_counts(grid.getDimX(), grid.getDimY(), grid.getDimZ(), 0);
        for (const auto& p : grid.particles()) {
            int i = static_cast<int>(p.position.x / grid.getDx());
            int j = static_cast<int>(p.position.y / grid.getDx());
            int k = static_cast<int>(p.position.z / grid.getDx());
            if (i >= 0 && i < grid.getDimX() && j >= 0 && j < grid.getDimY() && k >= 0 && k < grid.getDimZ()) {
                ++particle_counts(i, j, k);
            }
        }

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

        auto& particles = grid.particles();
        particles.erase(
            std::remove_if(
                particles.begin(),
                particles.end(),
                [&](const Particles& p) {
                    int i = static_cast<int>(p.position.x / grid.getDx());
                    int j = static_cast<int>(p.position.y / grid.getDx());
                    int k = static_cast<int>(p.position.z / grid.getDx());

                    if (i < 0 || i >= grid.getDimX() || j < 0 || j >= grid.getDimY() || k < 0 || k >= grid.getDimZ()) {
                        return true;
                    }

                    if (grid.celltypes()(i, j, k) == CellType::FLUID) {
                        int& c = particle_counts(i, j, k);
                        int& d = particles_to_drop(i, j, k);
                        if (d > 0 && c > 0) {
                            float prob = static_cast<float>(d) / static_cast<float>(c);
                            --c;
                            if (static_cast<float>(rand()) / RAND_MAX < prob) {
                                --d;
                                return true;
                            }
                        } else if (c > 0) {
                            --c;
                        }
                    }
                    return false;
                }),
            particles.end());

        mc.reset();
        for (int k = 0; k < grid.getDimZ(); ++k) {
            for (int j = 0; j < grid.getDimY(); ++j) {
                for (int i = 0; i < grid.getDimX(); ++i) {
                    mc.set_data(grid.liquid_phi()(i, j, k), i, j, k);
                }
            }
        }

        mc.set_method(false);
        mc.run(0.0f);

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

        std::stringstream filename;
        filename << "frame_" << std::setw(4) << std::setfill('0') << frame << ".obj";
        DataExporter::exportMeshToObj(render_mesh, (output_dir / filename.str()).string());
    }

    std::cout << "Simulation completed." << std::endl;
    return 0;
}
