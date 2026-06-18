#include "solver.h"

#include <algorithm>

namespace Solver {
    Float dotProduct(const Grid<Float>& a, const Grid<Float>& b, const Grid<CellType>& cellTypes) {
        Float result = 0.0f;
        for (int k = 0; k < a.getDepth(); ++k) {
            for (int j = 0; j < a.getHeight(); ++j) {
                for (int i = 0; i < a.getWidth(); ++i) {
                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        result += a(i, j, k) * b(i, j, k);
                    }
                }
            }
        }
        return result;
    }

    void saxpy(Grid<Float>& y, Float a, const Grid<Float>& x, const Grid<CellType>& cellTypes) {
        for (int k = 0; k < y.getDepth(); ++k) {
            for (int j = 0; j < y.getHeight(); ++j) {
                for (int i = 0; i < y.getWidth(); ++i) {
                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        y(i, j, k) += a * x(i, j, k);
                    }
                }
            }
        }
    }

    Grid<Float> discrete_divergence(const MACGrid& grid) {
        const int nx = grid.getDimX();
        const int ny = grid.getDimY();
        const int nz = grid.getDimZ();
        const Float dx = grid.getDx();
        Grid<Float> negative_divergence(nx, ny, nz, 0.0f);

        const auto& u = grid.u();
        const auto& v = grid.v();
        const auto& w = grid.w();
        const auto& celltypes = grid.celltypes();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (celltypes(i, j, k) == CellType::FLUID) {
                        Float divergence =
                            (u(i + 1, j, k) - u(i, j, k)) +
                            (v(i, j + 1, k) - v(i, j, k)) +
                            (w(i, j, k + 1) - w(i, j, k));
                        negative_divergence(i, j, k) = -divergence / dx;
                    }
                }
            }
        }

        return negative_divergence;
    }

    void buildMatrixA(
        Grid<Float>& Adiag,
        Grid<Float>& Aplus_i,
        Grid<Float>& Aplus_j,
        Grid<Float>& Aplus_k,
        const Grid<CellType>& cellTypes,
        Float dx,
        Float dt,
        Float rho) {
        const int nx = Adiag.getWidth();
        const int ny = Adiag.getHeight();
        const int nz = Adiag.getDepth();
        const Float scale = dt / (rho * dx * dx);

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;

                    if (cellTypes(i, j, k) != CellType::FLUID) {
                        continue;
                    }

                    Float diag_val = 0.0f;
                    if (i < nx - 1 && cellTypes(i + 1, j, k) != CellType::SOLID) ++diag_val;
                    if (i > 0 && cellTypes(i - 1, j, k) != CellType::SOLID) ++diag_val;
                    if (j < ny - 1 && cellTypes(i, j + 1, k) != CellType::SOLID) ++diag_val;
                    if (j > 0 && cellTypes(i, j - 1, k) != CellType::SOLID) ++diag_val;
                    if (k < nz - 1 && cellTypes(i, j, k + 1) != CellType::SOLID) ++diag_val;
                    if (k > 0 && cellTypes(i, j, k - 1) != CellType::SOLID) ++diag_val;

                    Adiag(i, j, k) = (diag_val == 0.0f) ? scale : diag_val * scale;
                    if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) Aplus_i(i, j, k) = -scale;
                    if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) Aplus_j(i, j, k) = -scale;
                    if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) Aplus_k(i, j, k) = -scale;
                }
            }
        }
    }

    void applyA(
        Grid<Float>& result,
        const Grid<Float>& p,
        const Grid<CellType>& cellTypes,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k) {
        const int nx = result.getWidth();
        const int ny = result.getHeight();
        const int nz = result.getDepth();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (cellTypes(i, j, k) != CellType::FLUID) {
                        result(i, j, k) = 0.0f;
                        continue;
                    }

                    Float val = Adiag(i, j, k) * p(i, j, k);
                    if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) val += Aplus_i(i, j, k) * p(i + 1, j, k);
                    if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) val += Aplus_i(i - 1, j, k) * p(i - 1, j, k);
                    if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) val += Aplus_j(i, j, k) * p(i, j + 1, k);
                    if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) val += Aplus_j(i, j - 1, k) * p(i, j - 1, k);
                    if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) val += Aplus_k(i, j, k) * p(i, j, k + 1);
                    if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) val += Aplus_k(i, j, k - 1) * p(i, j, k - 1);
                    result(i, j, k) = val;
                }
            }
        }
    }

    Float dotProduct_FVM(const Grid<Float>& a, const Grid<Float>& b, const MACGrid& grid) {
        Float result = 0.0f;
        const auto& celltypes = grid.celltypes();
        for (int k = 0; k < a.getDepth(); ++k) {
            for (int j = 0; j < a.getHeight(); ++j) {
                for (int i = 0; i < a.getWidth(); ++i) {
                    if (celltypes(i, j, k) == CellType::FLUID) {
                        result += a(i, j, k) * b(i, j, k);
                    }
                }
            }
        }
        return result;
    }

    void saxpy_FVM(Grid<Float>& y, Float a, const Grid<Float>& x, const MACGrid& grid) {
        const auto& celltypes = grid.celltypes();
        for (int k = 0; k < y.getDepth(); ++k) {
            for (int j = 0; j < y.getHeight(); ++j) {
                for (int i = 0; i < y.getWidth(); ++i) {
                    if (celltypes(i, j, k) == CellType::FLUID) {
                        y(i, j, k) += a * x(i, j, k);
                    }
                }
            }
        }
    }

    Grid<Float> discrete_divergence_FVM(const MACGrid& grid) {
        const int nx = grid.getDimX();
        const int ny = grid.getDimY();
        const int nz = grid.getDimZ();
        const Float dx = grid.getDx();
        Grid<Float> negative_divergence(nx, ny, nz, 0.0f);

        const auto& u = grid.u();
        const auto& v = grid.v();
        const auto& w = grid.w();
        const auto& u_area = grid.area_u();
        const auto& v_area = grid.area_v();
        const auto& w_area = grid.area_w();
        const auto& celltypes = grid.celltypes();
        const auto& solidVelocity = grid.solidvelocity();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (celltypes(i, j, k) != CellType::FLUID) {
                        continue;
                    }

                    Float fluid_flux =
                        (u(i + 1, j, k) * u_area(i + 1, j, k) - u(i, j, k) * u_area(i, j, k)) +
                        (v(i, j + 1, k) * v_area(i, j + 1, k) - v(i, j, k) * v_area(i, j, k)) +
                        (w(i, j, k + 1) * w_area(i, j, k + 1) - w(i, j, k) * w_area(i, j, k));

                    int i_p1 = std::min(i + 1, nx - 1);
                    int i_m1 = std::max(i - 1, 0);
                    int j_p1 = std::min(j + 1, ny - 1);
                    int j_m1 = std::max(j - 1, 0);
                    int k_p1 = std::min(k + 1, nz - 1);
                    int k_m1 = std::max(k - 1, 0);

                    Vector3f solid_vel_right = 0.5f * (solidVelocity(i_p1, j, k) + solidVelocity(i, j, k));
                    Vector3f solid_vel_left = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i_m1, j, k));
                    Vector3f solid_vel_top = 0.5f * (solidVelocity(i, j_p1, k) + solidVelocity(i, j, k));
                    Vector3f solid_vel_bottom = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i, j_m1, k));
                    Vector3f solid_vel_front = 0.5f * (solidVelocity(i, j, k_p1) + solidVelocity(i, j, k));
                    Vector3f solid_vel_back = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i, j, k_m1));

                    Float solid_flux = 0.0f;
                    solid_flux += solid_vel_right.x * (1.0f - u_area(i + 1, j, k));
                    solid_flux -= solid_vel_left.x * (1.0f - u_area(i, j, k));
                    solid_flux += solid_vel_top.y * (1.0f - v_area(i, j + 1, k));
                    solid_flux -= solid_vel_bottom.y * (1.0f - v_area(i, j, k));
                    solid_flux += solid_vel_front.z * (1.0f - w_area(i, j, k + 1));
                    solid_flux -= solid_vel_back.z * (1.0f - w_area(i, j, k));

                    negative_divergence(i, j, k) = -(fluid_flux + solid_flux) / dx;
                }
            }
        }

        return negative_divergence;
    }

    void buildMatrixA_FVM(
        Grid<Float>& Adiag,
        Grid<Float>& Aplus_i,
        Grid<Float>& Aplus_j,
        Grid<Float>& Aplus_k,
        const MACGrid& grid,
        const Grid<CellType>& cellTypes,
        Float dt,
        Float rho) {
        const int nx = grid.getDimX();
        const int ny = grid.getDimY();
        const int nz = grid.getDimZ();
        const Float dx = grid.getDx();
        const auto& u_area = grid.area_u();
        const auto& v_area = grid.area_v();
        const auto& w_area = grid.area_w();
        const auto& density = grid.density();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;

                    if (cellTypes(i, j, k) != CellType::FLUID) {
                        continue;
                    }

                    Float diag_val = 0.0f;

                    auto safe_rho = [&](Float r) { return (r < 1e-6f) ? rho : r; };

                    if (i < nx - 1) {
                        Float rho_face = safe_rho(0.5f * (density(i, j, k) + density(i + 1, j, k)));
                        Float term = dt / (rho_face * dx * dx) * u_area(i + 1, j, k);
                        if (cellTypes(i + 1, j, k) == CellType::FLUID) {
                            diag_val += term;
                            Aplus_i(i, j, k) = -term;
                        } else if (cellTypes(i + 1, j, k) != CellType::SOLID) {
                            diag_val += term;
                        }
                    }
                    if (i > 0) {
                        Float rho_face = safe_rho(0.5f * (density(i, j, k) + density(i - 1, j, k)));
                        Float term = dt / (rho_face * dx * dx) * u_area(i, j, k);
                        if (cellTypes(i - 1, j, k) != CellType::SOLID) {
                            diag_val += term;
                        }
                    }
                    if (j < ny - 1) {
                        Float rho_face = safe_rho(0.5f * (density(i, j, k) + density(i, j + 1, k)));
                        Float term = dt / (rho_face * dx * dx) * v_area(i, j + 1, k);
                        if (cellTypes(i, j + 1, k) == CellType::FLUID) {
                            diag_val += term;
                            Aplus_j(i, j, k) = -term;
                        } else if (cellTypes(i, j + 1, k) != CellType::SOLID) {
                            diag_val += term;
                        }
                    }
                    if (j > 0) {
                        Float rho_face = safe_rho(0.5f * (density(i, j, k) + density(i, j - 1, k)));
                        Float term = dt / (rho_face * dx * dx) * v_area(i, j, k);
                        if (cellTypes(i, j - 1, k) != CellType::SOLID) {
                            diag_val += term;
                        }
                    }
                    if (k < nz - 1) {
                        Float rho_face = safe_rho(0.5f * (density(i, j, k) + density(i, j, k + 1)));
                        Float term = dt / (rho_face * dx * dx) * w_area(i, j, k + 1);
                        if (cellTypes(i, j, k + 1) == CellType::FLUID) {
                            diag_val += term;
                            Aplus_k(i, j, k) = -term;
                        } else if (cellTypes(i, j, k + 1) != CellType::SOLID) {
                            diag_val += term;
                        }
                    }
                    if (k > 0) {
                        Float rho_face = safe_rho(0.5f * (density(i, j, k) + density(i, j, k - 1)));
                        Float term = dt / (rho_face * dx * dx) * w_area(i, j, k);
                        if (cellTypes(i, j, k - 1) != CellType::SOLID) {
                            diag_val += term;
                        }
                    }

                    Float safe_density = safe_rho(density(i, j, k));
                    Adiag(i, j, k) = (diag_val == 0.0f) ? dt / (safe_density * dx * dx) : diag_val;
                }
            }
        }
    }

    void applyA_FVM(
        Grid<Float>& result,
        const Grid<Float>& p,
        const MACGrid& grid,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k) {
        applyA(result, p, grid.celltypes(), Adiag, Aplus_i, Aplus_j, Aplus_k);
    }

    void MIC0preconditioner(
        Grid<Float>& precon,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const Grid<CellType>& cellTypes) {
        const int nx = precon.getWidth();
        const int ny = precon.getHeight();
        const int nz = precon.getDepth();
        const Float tau = 0.97f;
        const Float sgm = 0.25f;

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (cellTypes(i, j, k) != CellType::FLUID) {
                        precon(i, j, k) = 0.0f;
                        continue;
                    }

                    Float term_i_sq = 0.0f;
                    Float term_j_sq = 0.0f;
                    Float term_k_sq = 0.0f;
                    if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                        term_i_sq = Sqr(Aplus_i(i - 1, j, k) * precon(i - 1, j, k));
                    }
                    if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                        term_j_sq = Sqr(Aplus_j(i, j - 1, k) * precon(i, j - 1, k));
                    }
                    if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                        term_k_sq = Sqr(Aplus_k(i, j, k - 1) * precon(i, j, k - 1));
                    }

                    Float comp_i = 0.0f;
                    Float comp_j = 0.0f;
                    Float comp_k = 0.0f;
                    if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                        Float precon_sq = Sqr(precon(i - 1, j, k));
                        comp_i = Aplus_i(i - 1, j, k) * (Aplus_j(i - 1, j, k) + Aplus_k(i - 1, j, k)) * precon_sq;
                    }
                    if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                        Float precon_sq = Sqr(precon(i, j - 1, k));
                        comp_j = Aplus_j(i, j - 1, k) * (Aplus_i(i, j - 1, k) + Aplus_k(i, j - 1, k)) * precon_sq;
                    }
                    if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                        Float precon_sq = Sqr(precon(i, j, k - 1));
                        comp_k = Aplus_k(i, j, k - 1) * (Aplus_i(i, j, k - 1) + Aplus_j(i, j, k - 1)) * precon_sq;
                    }

                    Float e = Adiag(i, j, k) - term_i_sq - term_j_sq - term_k_sq - tau * (comp_i + comp_j + comp_k);
                    if (e < sgm * Adiag(i, j, k)) {
                        e = Adiag(i, j, k);
                    }
                    precon(i, j, k) = 1.0f / std::sqrt(e);
                }
            }
        }
    }

    void applyPreconditioner(
        Grid<Float>& z,
        const Grid<Float>& r,
        const Grid<Float>& precon,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const Grid<CellType>& cellTypes) {
        const int nx = z.getWidth();
        const int ny = z.getHeight();
        const int nz = z.getDepth();
        Grid<Float> q(nx, ny, nz, 0.0f);

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (cellTypes(i, j, k) != CellType::FLUID) {
                        continue;
                    }

                    Float offdiag_sum = 0.0f;
                    if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) offdiag_sum += Aplus_i(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
                    if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) offdiag_sum += Aplus_j(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
                    if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) offdiag_sum += Aplus_k(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                    q(i, j, k) = (r(i, j, k) - offdiag_sum) * precon(i, j, k);
                }
            }
        }

        for (int k = nz - 1; k >= 0; --k) {
            for (int j = ny - 1; j >= 0; --j) {
                for (int i = nx - 1; i >= 0; --i) {
                    if (cellTypes(i, j, k) != CellType::FLUID) {
                        continue;
                    }

                    Float offdiag_sum = 0.0f;
                    if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) offdiag_sum += Aplus_i(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
                    if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) offdiag_sum += Aplus_j(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
                    if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) offdiag_sum += Aplus_k(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                    z(i, j, k) = (q(i, j, k) - offdiag_sum) * precon(i, j, k);
                }
            }
        }
    }

    void MIC0preconditioner_FVM(
        Grid<Float>& precon,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const MACGrid& grid) {
        MIC0preconditioner(precon, Adiag, Aplus_i, Aplus_j, Aplus_k, grid.celltypes());
    }

    void applyPreconditioner_FVM(
        Grid<Float>& z,
        const Grid<Float>& r,
        const Grid<Float>& precon,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const MACGrid& grid) {
        applyPreconditioner(z, r, precon, Aplus_i, Aplus_j, Aplus_k, grid.celltypes());
    }

    void PCG(
        Grid<Float>& p,
        const Grid<Float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        Float dx,
        int maxIterations,
        Float tolerance) {
        (void)dx;
        const int nx = p.getWidth();
        const int ny = p.getHeight();
        const int nz = p.getDepth();

        Grid<Float> r(nx, ny, nz, 0.0f);
        Grid<Float> Ap(nx, ny, nz, 0.0f);
        r = b;
        applyA(Ap, p, cellTypes, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
        saxpy(r, -1.0f, Ap, cellTypes);

        Float initial_residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
        if (initial_residual_norm < 1e-9f) {
            std::cout << "Initial residual is zero." << std::endl;
            return;
        }

        Grid<Float> z(nx, ny, nz, 0.0f);
        applyPreconditioner(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
        Grid<Float> d = z;
        Float delta_new = dotProduct(r, z, cellTypes);

        for (int iter = 0; iter < maxIterations; ++iter) {
            Grid<Float> q(nx, ny, nz, 0.0f);
            applyA(q, d, cellTypes, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);

            Float d_q = dotProduct(d, q, cellTypes);
            if (std::abs(d_q) < static_cast<Float>(1e-20)) {
                std::cout << "PCG stopped: search direction is orthogonal to A*d." << std::endl;
                return;
            }

            Float alpha = delta_new / d_q;
            saxpy(p, alpha, d, cellTypes);
            saxpy(r, -alpha, q, cellTypes);

            Float residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
            std::cout << "Iteration " << iter + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "Converged after " << iter + 1 << " iterations." << std::endl;
                return;
            }

            applyPreconditioner(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
            Float delta_old = delta_new;
            delta_new = dotProduct(r, z, cellTypes);
            if (std::abs(delta_old) < static_cast<Float>(1e-20)) {
                std::cout << "PCG stopped: delta is too small." << std::endl;
                return;
            }

            Float beta = delta_new / delta_old;
            Grid<Float> temp_d = d;
            d = z;
            saxpy(d, beta, temp_d, cellTypes);
        }

        std::cout << "PCG did not converge after " << maxIterations << " iterations." << std::endl;
    }

    void PCG_FVM(
        Grid<Float>& p,
        const Grid<Float>& b,
        const SystemMatrix& matrix,
        const MACGrid& grid,
        int maxIterations,
        Float tolerance) {
        const int nx = p.getWidth();
        const int ny = p.getHeight();
        const int nz = p.getDepth();

        Grid<Float> r(nx, ny, nz, 0.0f);
        Grid<Float> Ap(nx, ny, nz, 0.0f);
        r = b;
        applyA_FVM(Ap, p, grid, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
        saxpy_FVM(r, -1.0f, Ap, grid);

        Float initial_residual_norm = std::sqrt(dotProduct_FVM(r, r, grid));
        if (initial_residual_norm < 1e-9f) {
            std::cout << "FVM: Initial residual is zero." << std::endl;
            return;
        }

        Grid<Float> z(nx, ny, nz, 0.0f);
        applyPreconditioner_FVM(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
        Grid<Float> d = z;
        Float delta_new = dotProduct_FVM(r, z, grid);

        for (int iter = 0; iter < maxIterations; ++iter) {
            Grid<Float> q(nx, ny, nz, 0.0f);
            applyA_FVM(q, d, grid, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);

            Float d_q = dotProduct_FVM(d, q, grid);
            if (std::abs(d_q) < static_cast<Float>(1e-20)) {
                std::cout << "FVM PCG stopped: search direction is orthogonal to A*d." << std::endl;
                return;
            }

            Float alpha = delta_new / d_q;
            saxpy_FVM(p, alpha, d, grid);
            saxpy_FVM(r, -alpha, q, grid);

            Float residual_norm = std::sqrt(dotProduct_FVM(r, r, grid));
            std::cout << "FVM Iteration " << iter + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "FVM Converged after " << iter + 1 << " iterations." << std::endl;
                return;
            }

            applyPreconditioner_FVM(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
            Float delta_old = delta_new;
            delta_new = dotProduct_FVM(r, z, grid);
            if (std::abs(delta_old) < static_cast<Float>(1e-20)) {
                std::cout << "FVM PCG stopped: delta is too small." << std::endl;
                return;
            }

            Float beta = delta_new / delta_old;
            Grid<Float> temp_d = d;
            d = z;
            saxpy_FVM(d, beta, temp_d, grid);
        }

        std::cout << "FVM PCG did not converge after " << maxIterations << " iterations." << std::endl;
    }
}
