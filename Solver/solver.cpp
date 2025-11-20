#include"solver.h"

namespace Solver {
    Float Solver::dotProduct(const Grid<Float>& a, const Grid<Float>& b, const Grid<CellType>& cellTypes) {
    Float result = 0.0f;
    int nx = a.getWidth();
    int ny = a.getHeight();
    int nz = a.getDepth();
    for(int k=0;k<nz;++k){
        for(int j=0;j<ny;++j){
            for(int i=0;i<nx;++i){
                if(cellTypes(i, j, k) == CellType::FLUID) {
                    result += a(i, j, k) * b(i, j, k);
                }
            }
        }
    }
    return result;
    }
    void Solver::saxpy(Grid<Float>& y, Float a, const Grid<Float>& x, const Grid<CellType>& cellTypes) {
        int nx = y.getWidth();
        int ny = y.getHeight();
        int nz = y.getDepth();
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i, j, k) == CellType::FLUID) {
                        y(i, j, k) += a * x(i, j, k);
                    }
                }
            }
        }
    }
    Grid<Float> Solver::discrete_divergence(const MACGrid& grid){
        int nx = grid.celltypes().getWidth();
        int ny = grid.celltypes().getHeight();
        int nz = grid.celltypes().getDepth();
        Float dx = grid.getDx();
        Grid<Float> negetivedivergence(nx,ny,nz,0.0f);
        const auto& u = grid.u();
        const auto& v = grid.v();
        const auto& w = grid.w();
        const auto& celltypes = grid.celltypes();
        for(int k = 0; k<nz;++k){
            for(int j = 0; j<ny;++j){
                for(int i = 0;i<nx;++i){
                    if(celltypes(i,j,k)==CellType::FLUID){
                        Float u_right = u(i+1,j,k);
                        Float u_left = u(i,j,k);
                        Float v_top = v(i,j+1,k);
                        Float v_bottom = v(i,j,k);
                        Float w_front = w(i,j,k+1);
                        Float w_back = w(i,j,k);
                        Float divergence = (u_right-u_left)+(v_top-v_bottom)+(w_front-w_back);
                        negetivedivergence(i,j,k)=-divergence/dx;
                    }
                }
            }
        }
        return negetivedivergence;
    }
    void Solver::buildMatrixA(
    Grid<Float>& Adiag, 
    Grid<Float>& Aplus_i, 
    Grid<Float>& Aplus_j, 
    Grid<Float>& Aplus_k,
    const Grid<CellType>& cellTypes,
    Float dx, 
    Float dt, 
    Float rho
    ) {
        int nx = Adiag.getWidth();
        int ny = Adiag.getHeight();
        int nz = Adiag.getDepth();
        Float scale = dt / (rho * dx * dx);

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    // 初始化，确保每次循环都是干净的
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;

                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        Float diag_val = 0;
                        if (i < nx - 1 && cellTypes(i + 1, j, k) != CellType::AIR) diag_val++;
                        if (i > 0     && cellTypes(i - 1, j, k) != CellType::AIR) diag_val++;
                        if (j < ny - 1 && cellTypes(i, j + 1, k) != CellType::AIR) diag_val++;
                        if (j > 0     && cellTypes(i, j - 1, k) != CellType::AIR) diag_val++;
                        if (k < nz - 1 && cellTypes(i, j, k + 1) != CellType::AIR) diag_val++;
                        if (k > 0     && cellTypes(i, j, k - 1) != CellType::AIR) diag_val++;
                        
                        Adiag(i,j,k) = diag_val * scale;

                        if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) Aplus_i(i,j,k) = -scale;
                        if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) Aplus_j(i,j,k) = -scale;
                        if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) Aplus_k(i,j,k) = -scale;
                    }
                }
            }
        }
    }
    void Solver::applyA(
        Grid<Float>& result, 
        const Grid<Float>& p,
        const Grid<CellType>& cellTypes,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k
    ) {
        int nx = result.getWidth();
        int ny = result.getHeight();
        int nz = result.getDepth();

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        
                        // result_i = Σ (A_ij * p_j)
                        //         = A_ii*p_i + Σ (A_in * p_n) for neighbors n

                        // 1. 对角线部分 A(i,j,k),(i,j,k) * p(i,j,k)
                        Float val = Adiag(i, j, k) * p(i, j, k);

                        // 2. 非对角线部分（邻居的影响）
                        // X方向
                        if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) {
                            val += Aplus_i(i, j, k) * p(i + 1, j, k);
                        }
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            // 利用对称性: A(i,i-1) = A(i-1,i) = Aplus_i(i-1)
                            val += Aplus_i(i - 1, j, k) * p(i - 1, j, k);
                        }

                        // Y方向
                        if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) {
                            val += Aplus_j(i, j, k) * p(i, j + 1, k);
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            val += Aplus_j(i, j - 1, k) * p(i, j - 1, k);
                        }

                        // Z方向
                        if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) {
                            val += Aplus_k(i, j, k) * p(i, j, k + 1);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            val += Aplus_k(i, j, k - 1) * p(i, j, k - 1);
                        }

                        result(i, j, k) = val;

                    } else {
                        result(i, j, k) = 0.0f; // 非流体单元格结果为零
                    }
                }
            }
        }
    }
    Float Solver::dotProduct_FVM(const Grid<Float>& a, const Grid<Float>& b, const MACGrid& grid) {
        Float result = 0.0f;
        int nx = a.getWidth();
        int ny = a.getHeight();
        int nz = a.getDepth();
        const auto& volumeFractions = grid.volumeFractions();
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (volumeFractions(i, j, k) > 0.0f) {
                        result += a(i, j, k) * b(i, j, k) * volumeFractions(i, j, k);
                    }
                }
            }
        }
        return result;
    }
    void Solver::saxpy_FVM(Grid<Float>& y, Float a, const Grid<Float>& x, const MACGrid& grid){
        int nx = y.getWidth();
        int ny = y.getHeight();
        int nz = y.getDepth();
        const auto& volumeFractions = grid.volumeFractions();
        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    //判断标准改为体积分数
                    if(volumeFractions(i, j, k) > 0.0f) {
                        y(i, j, k) += a * x(i, j, k);
                    }
                }
            }
        }
    }
    Grid<Float> Solver::discrete_divergence_FVM(const MACGrid& grid){
        int nx = grid.getDimX();
        int ny = grid.getDimY();
        int nz = grid.getDimZ();
        Float dx = grid.getDx();
        Grid<Float> negetivedivergence(nx, ny, nz, 0.0f);
        
        const auto& u = grid.u();
        const auto& v = grid.v();
        const auto& w = grid.w();
        
        const auto& u_area = grid.area_u();
        const auto& v_area = grid.area_v();
        const auto& w_area = grid.area_w();
        
        const auto& volumeFractions = grid.volumeFractions();
        const auto& solidVelocity = grid.solidvelocity();
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (volumeFractions(i, j, k) > 0.0f) {
                        Float fluid_flux = 
                            (u(i + 1, j, k) * u_area(i + 1, j, k) - u(i, j, k) * u_area(i, j, k)) +
                            (v(i, j + 1, k) * v_area(i, j + 1, k) - v(i, j, k) * v_area(i, j, k)) +
                            (w(i, j, k + 1) * w_area(i, j, k + 1) - w(i, j, k) * w_area(i, j, k));
                        Float solid_flux = 0.0f;
                        // X方向
                        //进行边界检查
                        int i_p1 = std::min(i + 1, nx - 1);
                        int i_m1 = std::max(i - 1, 0);
                        Vector3f solid_vel_right = 0.5f * (solidVelocity(i_p1, j, k) + solidVelocity(i, j, k));
                        Vector3f solid_vel_left  = 0.5f * (solidVelocity(i,j,k) + solidVelocity(i_m1,j,k));
                        solid_flux += solid_vel_right.x * (1.0f - u_area(i + 1, j, k));
                        solid_flux -= solid_vel_left.x  * (1.0f - u_area(i, j, k));
                        // Y方向
                        int j_p1 = std::min(j + 1, ny - 1);
                        int j_m1 = std::max(j - 1, 0);
                        Vector3f solid_vel_top    = 0.5f * (solidVelocity(i, j_p1, k) + solidVelocity(i, j, k));
                        Vector3f solid_vel_bottom = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i, j_m1, k));
                        solid_flux += solid_vel_top.y    * (1.0f - v_area(i, j + 1, k));
                        solid_flux -= solid_vel_bottom.y * (1.0f - v_area(i, j, k));

                        // Z方向
                        int k_p1 = std::min(k + 1, nz - 1);
                        int k_m1 = std::max(k - 1, 0);
                        Vector3f solid_vel_front = 0.5f * (solidVelocity(i, j, k_p1) + solidVelocity(i, j, k));
                        Vector3f solid_vel_back = 0.5f * (solidVelocity(i, j, k) + solidVelocity(i, j, k_m1));
                        solid_flux += solid_vel_front.z * (1.0f - w_area(i, j, k + 1));
                        solid_flux -= solid_vel_back.z * (1.0f - w_area(i, j, k));
                        // 计算负散度
                        negetivedivergence(i, j, k) = 
                            -(fluid_flux - solid_flux) / dx;
                    }
                }
            }
        }
        return negetivedivergence;
    };
    void buildMatrixA_FVM(
        Grid<Float>& Adiag, 
        Grid<Float>& Aplus_i, 
        Grid<Float>& Aplus_j, 
        Grid<Float>& Aplus_k,
        const MACGrid& grid,
        Float dt
    ){
        int nx = grid.getDimX();
        int ny = grid.getDimY();
        int nz = grid.getDimZ();
        Float dx = grid.getDx();
        const auto& u_area = grid.area_u();
        const auto& v_area = grid.area_v();
        const auto& w_area = grid.area_w();
        const auto& volume_frac = grid.volumeFractions();
        const auto& density = grid.density();
        for (int k = 0; k < nz; ++k){
            for (int j = 0; j < ny; ++j){
                for (int i = 0; i < nx; ++i){
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;
                    if (volume_frac(i ,j, k) > 0.0f){
                        Float diag_val = 0.0f;
                        //x方向
                        // 右侧
                        if (i < nx -1){
                            Float rho_face = 0.5 * (density(i, j, k) + density(i + 1, j, k));
                            Float scale_face = dt / (rho_face * dx * dx);
                            Aplus_i(i, j, k) = -scale_face * u_area(i + 1, j, k);
                            diag_val += scale_face * u_area(i + 1, j, k);
                        }
                        // 左侧
                        if (i > 0){
                            Float rho_face = 0.5 * (density(i, j, k) + density(i - 1, j, k));
                            Float scale_face = dt / (rho_face * dx * dx);
                            diag_val += scale_face * u_area(i, j, k);
                        }
                        //y方向
                        // 上侧
                        if (j < ny - 1) {
                            Float rho_face = 0.5 * (density(i, j, k) + density(i, j + 1, k));
                            Float scale_face = dt / (rho_face * dx * dx);
                            Aplus_j(i, j, k) = -scale_face * v_area(i, j + 1, k);
                            diag_val += scale_face * v_area(i, j + 1, k);
                        }
                        // 下侧
                        if (j > 0) {
                            Float rho_face = 0.5 * (density(i, j, k) + density(i, j - 1, k));
                            Float scale_face = dt / (rho_face * dx * dx);
                            diag_val += scale_face * v_area(i, j, k);
                        }
                        //z方向
                        // 前侧
                        if (k < nz - 1) {
                            Float rho_face = 0.5 * (density(i, j, k) + density(i, j, k + 1));
                            Float scale_face = dt / (rho_face * dx * dx);
                            Aplus_k(i, j, k) = -scale_face * w_area(i, j, k + 1);
                            diag_val += scale_face * w_area(i, j, k + 1);
                        }
                        // 后侧
                        if (k > 0) {
                            Float rho_face = 0.5 * (density(i, j, k) + density(i, j, k - 1));
                            Float scale_face = dt / (rho_face * dx * dx);
                            diag_val += scale_face * w_area(i, j, k);
                        }
                        Adiag(i, j, k) = diag_val;
                    }
                }
            }
        };
    }
    void Solver::applyA_FVM(
        Grid<Float>& result, 
        const Grid<Float>& p,
        const MACGrid& grid,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k
    ) {
        int nx = result.getWidth();
        int ny = result.getHeight();
        int nz = result.getDepth();
        const auto& volumeFraction = grid.volumeFractions();
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    if (grid.volumeFractions()(i, j, k) > 0.0f) {
                        // 1. 对角线部分 A(i,j,k),(i,j,k) * p(i,j,k)
                        Float val = Adiag(i, j, k) * p(i, j, k);

                        // 2. 非对角线部分（邻居的影响）
                        // X方向
                        if (i < nx - 1 && volumeFraction(i + 1, j, k) > 0.0f) {
                        val += Aplus_i(i, j, k) * p(i + 1, j, k);
                        }
                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            val += Aplus_i(i - 1, j, k) * p(i - 1, j, k);
                        }

                        // Y方向
                        if (j < ny - 1 && volumeFraction(i, j + 1, k) > 0.0f) {
                            val += Aplus_j(i, j, k) * p(i, j + 1, k);
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            val += Aplus_j(i, j - 1, k) * p(i, j - 1, k);
                        }

                        // Z方向
                        if (k < nz - 1 && volumeFraction(i, j, k + 1) > 0.0f) {
                            val += Aplus_k(i, j, k) * p(i, j, k + 1);
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            val += Aplus_k(i, j, k - 1) * p(i, j, k - 1);
                        }

                        result(i, j, k) = val;

                    } else {
                        result(i, j, k) = 0.0f; // 非流体单元格结果为零
                    }
                }
            }
        }
    }
    void MIC0preconditioner(
        Grid<Float>& precon,
        const Grid<Float>& Adiag,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const Grid<CellType>& cellTypes
    ){
        int nx = precon.getWidth();
        int ny = precon.getHeight();
        int nz = precon.getDepth();
        const Float tau = 0.97f; // MIC(0)的调节参数
        const Float sgm = 0.25f; // 安全参数
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        //IC(0)
                        Float term_i_sq = 0.0f, term_j_sq = 0.0f, term_k_sq = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            Float precon_val = precon(i - 1, j, k);
                            term_i_sq = (Aplus_i(i - 1, j, k) * precon_val) * (Aplus_i(i - 1, j, k) * precon_val);
                        }
                        if (j> 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            Float precon_val = precon(i, j - 1, k);
                            term_j_sq = (Aplus_j(i, j - 1, k) * precon_val) * (Aplus_j(i, j - 1, k) * precon_val);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            Float precon_val = precon(i, j, k - 1);
                            term_k_sq = (Aplus_k(i, j, k - 1) * precon_val) * (Aplus_k(i, j, k - 1) * precon_val);
                        }
                        //MIC(0)
                        Float comp_i = 0.0f, comp_j = 0.0f, comp_k = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            Float precon_sq = precon(i - 1, j, k) * precon(i - 1, j, k);
                            comp_i = Aplus_i(i - 1, j, k) * (Aplus_j(i - 1, j, k) + Aplus_k(i - 1, j, k)) * precon_sq;
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            Float precon_sq = precon(i, j - 1, k) * precon(i, j - 1, k);
                            comp_j = Aplus_j(i, j - 1, k) * (Aplus_i(i, j - 1, k) + Aplus_k(i, j - 1, k)) * precon_sq;
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            Float precon_sq = precon(i, j, k - 1) * precon(i, j, k - 1);
                            comp_k = Aplus_k(i, j, k - 1) * (Aplus_i(i, j, k - 1) + Aplus_j(i, j, k - 1)) * precon_sq;
                        }
                        Float e = Adiag(i, j, k) - term_i_sq - term_j_sq - term_k_sq - tau * (comp_i + comp_j + comp_k);
                        if (e < sgm * Adiag(i, j, k)) {
                        e = Adiag(i, j, k);
                        }
                        // 存储e的平方根的倒数
                        precon(i, j, k) = 1.0f / std::sqrt(e);
                    } else {
                        precon(i, j, k) = 0.0f;
                    }
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
        const Grid<CellType>& cellTypes
    ){
        int nx = z.getWidth();
        int ny = z.getHeight();
        int nz = z.getDepth();
        Grid<Float> q(nx, ny, nz, 0.0f);
        //前向替换 (Forward Substitution), 求解 Lq = r
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        Float offdiag_sum = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_i(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_j(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            offdiag_sum += Aplus_k(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                        }
                        Float t = r(i, j, k) - offdiag_sum;
                        q(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }
        //后向替换 (Backward Substitution), 求解 L^Tz = q
        for(int k=nz-1;k>=0;--k){
            for(int j=ny-1;j>=0;--j){
                for(int i=nx-1;i>=0;--i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        Float offdiag_sum = 0.0f;
                        if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_i(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
                        }
                        if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_j(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
                        }
                        if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) {
                            offdiag_sum += Aplus_k(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                        }
                        Float t = q(i, j, k) - offdiag_sum;
                        z(i, j, k) = t * precon(i, j, k);
                    }
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
        const MACGrid& grid // 传入 MACGrid 以获取 volumeFractions
    ){
        int nx = precon.getWidth();
        int ny = precon.getHeight();
        int nz = precon.getDepth();
        const Float tau = 0.97f;
        const Float sgm = 0.25f;

        const auto& volumeFraction = grid.volumeFractions();

        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    //判断标准改为体积分数
                    if(volumeFraction(i, j, k) > 0.0f){
                        
                        Float term_i_sq = 0.0f, term_j_sq = 0.0f, term_k_sq = 0.0f;

                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            Float precon_val = precon(i - 1, j, k);
                            term_i_sq = (Aplus_i(i - 1, j, k) * precon_val) * (Aplus_i(i - 1, j, k) * precon_val);
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            Float precon_val = precon(i, j - 1, k);
                            term_j_sq = (Aplus_j(i, j - 1, k) * precon_val) * (Aplus_j(i, j - 1, k) * precon_val);
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            Float precon_val = precon(i, j, k - 1);
                            term_k_sq = (Aplus_k(i, j, k - 1) * precon_val) * (Aplus_k(i, j, k - 1) * precon_val);
                        }
                        
                        Float comp_i = 0.0f, comp_j = 0.0f, comp_k = 0.0f;
                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            Float precon_sq = precon(i - 1, j, k) * precon(i - 1, j, k);
                            comp_i = Aplus_i(i - 1, j, k) * (Aplus_j(i - 1, j, k) + Aplus_k(i - 1, j, k)) * precon_sq;
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            Float precon_sq = precon(i, j - 1, k) * precon(i, j - 1, k);
                            comp_j = Aplus_j(i, j - 1, k) * (Aplus_i(i, j - 1, k) + Aplus_k(i, j - 1, k)) * precon_sq;
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            Float precon_sq = precon(i, j, k - 1) * precon(i, j, k - 1);
                            comp_k = Aplus_k(i, j, k - 1) * (Aplus_i(i, j, k - 1) + Aplus_j(i, j, k - 1)) * precon_sq;
                        }

                        Float e = Adiag(i, j, k) - term_i_sq - term_j_sq - term_k_sq - tau * (comp_i + comp_j + comp_k);
                        
                        if (e < sgm * Adiag(i, j, k)) {
                            e = Adiag(i, j, k);
                        }
                        
                        precon(i, j, k) = 1.0f / std::sqrt(e);

                    } else {
                        precon(i, j, k) = 0.0f;
                    }
                }
            }
        }
    }
    void Solver::applyPreconditioner_FVM(
        Grid<Float>& z,
        const Grid<Float>& r,
        const Grid<Float>& precon,
        const Grid<Float>& Aplus_i,
        const Grid<Float>& Aplus_j,
        const Grid<Float>& Aplus_k,
        const MACGrid& grid
    ){
        int nx = z.getWidth();
        int ny = z.getHeight();
        int nz = z.getDepth();
        Grid<Float> q(nx, ny, nz, 0.0f);
        
        const auto& volumeFraction = grid.volumeFractions();

        // 前向替换, 求解 Lq = r
        for(int k = 0; k < nz; ++k){
            for(int j = 0; j < ny; ++j){
                for(int i = 0; i < nx; ++i){
                    //判断标准改为体积分数
                    if(volumeFraction(i, j, k) > 0.0f){
                        Float offdiag_sum = 0.0f;
                        if (i > 0 && volumeFraction(i - 1, j, k) > 0.0f) {
                            offdiag_sum += Aplus_i(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
                        }
                        if (j > 0 && volumeFraction(i, j - 1, k) > 0.0f) {
                            offdiag_sum += Aplus_j(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
                        }
                        if (k > 0 && volumeFraction(i, j, k - 1) > 0.0f) {
                            offdiag_sum += Aplus_k(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                        }
                        Float t = r(i, j, k) - offdiag_sum;
                        q(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }

        // 后向替换, 求解 L^Tz = q
        for(int k = nz - 1; k >= 0; --k){
            for(int j = ny - 1; j >= 0; --j){
                for(int i = nx - 1; i >= 0; --i){
                    if(volumeFraction(i, j, k) > 0.0f){
                        Float offdiag_sum = 0.0f;
                        if (i < nx - 1 && volumeFraction(i + 1, j, k) > 0.0f) {
                            offdiag_sum += Aplus_i(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
                        }
                        if (j < ny - 1 && volumeFraction(i, j + 1, k) > 0.0f) {
                            offdiag_sum += Aplus_j(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
                        }
                        if (k < nz - 1 && volumeFraction(i, j, k + 1) > 0.0f) {
                            offdiag_sum += Aplus_k(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                        }
                        Float t = q(i, j, k) - offdiag_sum;
                        z(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }
    }
    void PCG(
        Grid<Float>& p,
        const Grid<Float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        Float dx,
        int maxIterations,
        Float tolerance
    ){
        //获取网格尺寸
        int nx = p.getWidth();
        int ny = p.getHeight();
        int nz = p.getDepth();
        //初始化
        Grid<Float> r(nx, ny, nz, 0.0f);
        Grid<Float> Ap(nx, ny, nz, 0.0f);
        r = b;
        //计算初始残差
        Float initial_residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
        if (initial_residual_norm < 1e-9) {
            std::cout << "Initial residual is zero."<<std::endl;
            return;
        }
        //z=M-1 * r
        Grid<Float> z(nx, ny, nz, 0.0f);
        applyPreconditioner(z,r,matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
        //d = z
        Grid<Float> d = z;
        // delta_new = r · z
        Float delta_new = dotProduct(r, z, cellTypes);

        //PCG
        for(int k =0; k < maxIterations;++k){
            //q = A * d
            Grid<Float> q(nx, ny, nz, 0.0f);
            applyA(q, d, cellTypes, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
            //alpha = delta_new / (d · q)
            Float d_q = dotProduct(d, q, cellTypes);
            Float alpha = delta_new / d_q;
            //p = p + alpha * d
            saxpy(p, alpha, d, cellTypes);
            //r = r - alpha * q
            saxpy(r, -alpha, q, cellTypes);
            //收敛性检查
            Float residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
            std::cout << "Iteration " << k + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
                return;
            }
            //z_new = M-1 * r_new
            applyPreconditioner(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
            Float delta_old = delta_new;
            delta_new = dotProduct(r, z, cellTypes);
            //beta = delta_new / delta_old
            Float beta = delta_new / delta_old;
            //d = z + beta * d
            Grid<Float> temp_d = d;
            d = z;
            saxpy(d, beta, temp_d, cellTypes);
        }
        std::cout << "PCG did not converge after" << maxIterations << std::endl;
    };
    void Solver::PCG_FVM(
        Grid<Float>& p,
        const Grid<Float>& b,
        const SystemMatrix& matrix,
        const MACGrid& grid,
        int maxIterations,
        Float tolerance
    ){
        int nx = p.getWidth();
        int ny = p.getHeight();
        int nz = p.getDepth();
        // 初始化残差
        Grid<Float> r(nx, ny, nz, 0.0f);
        r = b;
        // 计算初始残差范数用于收敛判断
        //调用 FVM 版本的dotProduct
        Float initial_residual_norm = std::sqrt(dotProduct_FVM(r, r, grid));
        if (initial_residual_norm < 1e-9) {
            std::cout << "FVM: Initial residual is zero." << std::endl;
            return;
        }
        // z = M^-1 * r
        Grid<Float> z(nx, ny, nz, 0.0f);
        //调用 FVM 版本的applyPreconditioner
        applyPreconditioner_FVM(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
        // d = z
        Grid<Float> d = z;
        // delta_new = r · z
        Float delta_new = dotProduct_FVM(r, z, grid);
        // --- 主循环 ---
        for (int k = 0; k < maxIterations; ++k) {
            // q = A * d
            Grid<Float> q(nx, ny, nz, 0.0f);
            applyA_FVM(q, d, grid, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
            
            // alpha = delta_new / (d · q)
            Float d_q = dotProduct_FVM(d, q, grid);
            Float alpha = delta_new / d_q;

            // p = p + alpha * d
            saxpy_FVM(p, alpha, d, grid);

            // r = r - alpha * q
            saxpy_FVM(r, -alpha, q, grid);

            // 收敛性检查
            Float residual_norm = std::sqrt(dotProduct_FVM(r, r, grid));
            std::cout << "FVM Iteration " << k + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "FVM Converged after " << k + 1 << " iterations." << std::endl;
                return;
            }
            
            // z_new = M^-1 * r_new
            applyPreconditioner_FVM(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, grid);
            
            Float delta_old = delta_new;
            //调用 FVM 版本的 dotProduct
            delta_new = dotProduct_FVM(r, z, grid);

            // beta = delta_new / delta_old
            Float beta = delta_new / delta_old;
            
            // d = z + beta * d
            Grid<Float> temp_d = d;
            d = z;
            saxpy_FVM(d, beta, temp_d, grid);
        }
        
        std::cout << "FVM PCG did not converge after " << maxIterations << " iterations." << std::endl;
    }
}