#include"solver.h"
#include<vector>

namespace Solver {
    float dotProduct(const Grid<float>& a, const Grid<float>& b, const Grid<CellType>& cellTypes) {
    float result = 0.0f;
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
    void saxpy(Grid<float>& y, float a, const Grid<float>& x, const Grid<CellType>& cellTypes) {
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
    void buildMatrixA(
    Grid<float>& Adiag, 
    Grid<float>& Aplus_i, 
    Grid<float>& Aplus_j, 
    Grid<float>& Aplus_k,
    const Grid<CellType>& cellTypes,
    float dx, 
    float dt, 
    float rho
    ) {
        int nx = Adiag.getWidth();
        int ny = Adiag.getHeight();
        int nz = Adiag.getDepth();
        float scale = dt / (rho * dx * dx);

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    // 初始化，确保每次循环都是干净的
                    Adiag(i, j, k) = 0.0f;
                    Aplus_i(i, j, k) = 0.0f;
                    Aplus_j(i, j, k) = 0.0f;
                    Aplus_k(i, j, k) = 0.0f;

                    if (cellTypes(i, j, k) == CellType::FLUID) {
                        float diag_val = 0;
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
    void applyA(
    Grid<float>& result, 
    const Grid<float>& p,
    const Grid<CellType>& cellTypes,
    const Grid<float>& Adiag,
    const Grid<float>& Aplus_i,
    const Grid<float>& Aplus_j,
    const Grid<float>& Aplus_k
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
                        float val = Adiag(i, j, k) * p(i, j, k);

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
    void MIC0preconditioner(
        Grid<float>& precon,
        const Grid<float>& Adiag,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes
    ){
        int nx = precon.getWidth();
        int ny = precon.getHeight();
        int nz = precon.getDepth();
        const float tau = 0.97f; // MIC(0)的调节参数
        const float sgm = 0.25f; // 安全参数
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        //IC(0)
                        float term_i_sq = 0.0f, term_j_sq = 0.0f, term_k_sq = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            float precon_val = precon(i - 1, j, k);
                            term_i_sq = (Aplus_i(i - 1, j, k) * precon_val) * (Aplus_i(i - 1, j, k) * precon_val);
                        }
                        if (j> 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            float precon_val = precon(i, j - 1, k);
                            term_j_sq = (Aplus_j(i, j - 1, k) * precon_val) * (Aplus_j(i, j - 1, k) * precon_val);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            float precon_val = precon(i, j, k - 1);
                            term_k_sq = (Aplus_k(i, j, k - 1) * precon_val) * (Aplus_k(i, j, k - 1) * precon_val);
                        }
                        //MIC(0)
                        float comp_i = 0.0f, comp_j = 0.0f, comp_k = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            float precon_sq = precon(i - 1, j, k) * precon(i - 1, j, k);
                            comp_i = Aplus_i(i - 1, j, k) * (Aplus_j(i - 1, j, k) + Aplus_k(i - 1, j, k)) * precon_sq;
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            float precon_sq = precon(i, j - 1, k) * precon(i, j - 1, k);
                            comp_j = Aplus_j(i, j - 1, k) * (Aplus_i(i, j - 1, k) + Aplus_k(i, j - 1, k)) * precon_sq;
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            float precon_sq = precon(i, j, k - 1) * precon(i, j, k - 1);
                            comp_k = Aplus_k(i, j, k - 1) * (Aplus_i(i, j, k - 1) + Aplus_j(i, j, k - 1)) * precon_sq;
                        }
                        float e = Adiag(i, j, k) - term_i_sq - term_j_sq - term_k_sq - tau * (comp_i + comp_j + comp_k);
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
        Grid<float>& z,
        const Grid<float>& r,
        const Grid<float>& precon,
        const Grid<float>& Aplus_i,
        const Grid<float>& Aplus_j,
        const Grid<float>& Aplus_k,
        const Grid<CellType>& cellTypes
    ){
        int nx = z.getWidth();
        int ny = z.getHeight();
        int nz = z.getDepth();
        Grid<float> q(nx, ny, nz, 0.0f);
        //前向替换 (Forward Substitution), 求解 Lq = r
        for(int k=0;k<nz;++k){
            for(int j=0;j<ny;++j){
                for(int i=0;i<nx;++i){
                    if(cellTypes(i,j,k)==CellType::FLUID){
                        float offdiag_sum = 0.0f;
                        if (i > 0 && cellTypes(i - 1, j, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_i(i - 1, j, k) * precon(i - 1, j, k) * q(i - 1, j, k);
                        }
                        if (j > 0 && cellTypes(i, j - 1, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_j(i, j - 1, k) * precon(i, j - 1, k) * q(i, j - 1, k);
                        }
                        if (k > 0 && cellTypes(i, j, k - 1) == CellType::FLUID) {
                            offdiag_sum += Aplus_k(i, j, k - 1) * precon(i, j, k - 1) * q(i, j, k - 1);
                        }
                        float t = r(i, j, k) - offdiag_sum;
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
                        float offdiag_sum = 0.0f;
                        if (i < nx - 1 && cellTypes(i + 1, j, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_i(i, j, k) * precon(i, j, k) * z(i + 1, j, k);
                        }
                        if (j < ny - 1 && cellTypes(i, j + 1, k) == CellType::FLUID) {
                            offdiag_sum += Aplus_j(i, j, k) * precon(i, j, k) * z(i, j + 1, k);
                        }
                        if (k < nz - 1 && cellTypes(i, j, k + 1) == CellType::FLUID) {
                            offdiag_sum += Aplus_k(i, j, k) * precon(i, j, k) * z(i, j, k + 1);
                        }
                        float t = q(i, j, k) - offdiag_sum;
                        z(i, j, k) = t * precon(i, j, k);
                    }
                }
            }
        }
    }
    void PCG(
        Grid<float>& p,
        const Grid<float>& b,
        const SystemMatrix& matrix,
        const Grid<CellType>& cellTypes,
        float dx,
        int maxIterations,
        float tolerance
    ){
        //获取网格尺寸
        int nx = p.getWidth();
        int ny = p.getHeight();
        int nz = p.getDepth();
        //初始化
        Grid<float> r(nx, ny, nz, 0.0f);
        Grid<float> Ap(nx, ny, nz, 0.0f);
        r = b;
        //计算初始残差
        float initial_residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
        if (initial_residual_norm < 1e-9) {
            std::cout << "Initial residual is zero."<<std::endl;
            return;
        }
        //z=M-1 * r
        Grid<float> z(nx, ny, nz, 0.0f);
        applyPreconditioner(z,r,matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
        //d = z
        Grid<float> d = z;
        // delta_new = r · z
        float delta_new = dotProduct(r, z, cellTypes);

        //PCG
        for(int k =0; k < maxIterations;++k){
            //q = A * d
            Grid<float> q(nx, ny, nz, 0.0f);
            applyA(q, d, cellTypes, matrix.Adiag, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k);
            //alpha = delta_new / (d · q)
            float d_q = dotProduct(d, q, cellTypes);
            float alpha = delta_new / d_q;
            //p = p + alpha * d
            saxpy(p, alpha, d, cellTypes);
            //r = r - alpha * q
            saxpy(r, -alpha, q, cellTypes);
            //收敛性检查
            float residual_norm = std::sqrt(dotProduct(r, r, cellTypes));
            std::cout << "Iteration " << k + 1 << ": Residual norm = " << residual_norm << std::endl;
            if (residual_norm < tolerance * initial_residual_norm) {
                std::cout << "Converged after " << k + 1 << " iterations." << std::endl;
                return;
            }
            //z_new = M-1 * r_new
            applyPreconditioner(z, r, matrix.precon, matrix.Aplus_i, matrix.Aplus_j, matrix.Aplus_k, cellTypes);
            float delta_old = delta_new;
            delta_new = dotProduct(r, z, cellTypes);
            //beta = delta_new / delta_old
            float beta = delta_new / delta_old;
            //d = z + beta * d
            Grid<float> temp_d = d;
            d = z;
            saxpy(d, beta, temp_d, cellTypes);
        }
        std::cout << "PCG did not converge after" << maxIterations << std::endl;
    };
}