#include "FLIP.h"
#include "advection.h"
namespace FLIPSolver {
    void ParticleToGrid(MACGrid& grid){
        //1. 获取u, v, w的网格引用, 并初始化为0
        auto& u = grid.u();
        auto& v = grid.v();
        auto& w = grid.w();
        u.fill(0.0f);
        v.fill(0.0f);
        w.fill(0.0f);
        //创建三个权重网格, 并初始化为0
        Grid<float> u_weight(u.getWidth(), u.getHeight(), u.getDepth(), 0.0f);
        Grid<float> v_weight(v.getWidth(), v.getHeight(), v.getDepth(), 0.0f);
        Grid<float> w_weight(w.getWidth(), w.getHeight(), w.getDepth(), 0.0f);
        const float dx = grid.getDx();
        //2. 遍历所有粒子
        for (const auto& particle : grid.particles()) {
            //获取粒子位置u方向
            const Point3f& u_pos = particle.position / dx;
            //计算粒子在网格中的索引
            int ui = static_cast<int>(u_pos.x);
            int uj = static_cast<int>(u_pos.y-0.5f);
            int uk = static_cast<int>(u_pos.z-0.5f);//这是在干什么-0.5?
            for(int k = 0; k<= 1; ++k){
                for(int j = 0; j<= 1; ++j){
                    for(int i = 0; i<= 1; ++i){
                        int i_ = ui + i;
                        int j_ = uj + j;
                        int k_ = uk + k;
                        // 计算三线性插值权重(以后可改用B样条)
                        float w_x = (i == 0) ? 1.0f - (u_pos.x - ui) : (u_pos.x - ui);
                        float w_y = (j == 0) ? 1.0f - (u_pos.y - 0.5f - uj) : (u_pos.y - 0.5f - uj);
                        float w_z = (k == 0) ? 1.0f - (u_pos.z - 0.5f - uk) : (u_pos.z - 0.5f - uk);
                        float weight = w_x * w_y * w_z;
                        //安全检查
                        if(i_>= 0 && i_ <u.getWidth() && j_>= 0 && j_ < u.getHeight() && k_ >= 0 && k_ < u.getDepth()){
                            u(i_, j_, k_) += particle.velocity.x * weight;
                            u_weight(i_, j_, k_) += weight;
                        }
                    }
                }
            }
                // --- 对V速度分量进行贡献 ---
            Point3f v_pos = particle.position / dx;
            int vi = static_cast<int>(v_pos.x - 0.5f);
            int vj = static_cast<int>(v_pos.y);
            int vk = static_cast<int>(v_pos.z - 0.5f);

            for (int k = 0; k <= 1; ++k) {
                for (int j = 0; j <= 1; ++j) {
                    for (int i = 0; i <= 1; ++i) {
                        int i_ = vi + i;
                        int j_ = vj + j;
                        int k_ = vk + k;

                        float w_x = (i == 0) ? 1.0f - (v_pos.x - 0.5f - vi) : (v_pos.x - 0.5f - vi);
                        float w_y = (j == 0) ? 1.0f - (v_pos.y - vj) : (v_pos.y - vj);
                        float w_z = (k == 0) ? 1.0f - (v_pos.z - 0.5f - vk) : (v_pos.z - 0.5f - vk);
                        float weight = w_x * w_y * w_z;

                        if (i_ >= 0 && i_ < v.getWidth() && j_ >= 0 && j_ < v.getHeight() && k_ >= 0 && k_ < v.getDepth()) {
                            v(i_, j_, k_) += particle.velocity.y * weight;
                            v_weight(i_, j_, k_) += weight;
                        }
                    }
                }
            }
        
            // --- 对W速度分量进行贡献 ---
            Point3f w_pos = particle.position / dx;
            int wi = static_cast<int>(w_pos.x - 0.5f);
            int wj = static_cast<int>(w_pos.y - 0.5f);
            int wk = static_cast<int>(w_pos.z);
            
            for (int k = 0; k <= 1; ++k) {
                for (int j = 0; j <= 1; ++j) {
                    for (int i = 0; i <= 1; ++i) {
                        int i_ = wi + i;
                        int j_ = wj + j;
                        int k_ = wk + k;

                        float w_x = (i == 0) ? 1.0f - (w_pos.x - 0.5f - wi) : (w_pos.x - 0.5f - wi);
                        float w_y = (j == 0) ? 1.0f - (w_pos.y - 0.5f - wj) : (w_pos.y - 0.5f - wj);
                        float w_z = (k == 0) ? 1.0f - (w_pos.z - wk) : (w_pos.z - wk);
                        float weight = w_x * w_y * w_z;

                        if (i_ >= 0 && i_ < w.getWidth() && j_ >= 0 && j_ < w.getHeight() && k_ >= 0 && k_ < w.getDepth()) {
                            w(i_, j_, k_) += particle.velocity.z * weight;
                            w_weight(i_, j_, k_) += weight;
                        }
                    }
                }
            }
        }
        // 4. 归一化，得到加权平均速度
        for (int k = 0; k < u.getDepth(); ++k) for (int j = 0; j < u.getHeight(); ++j) for (int i = 0; i < u.getWidth(); ++i) {
            if (u_weight(i, j, k) > 1e-9) u(i, j, k) /= u_weight(i, j, k);
        }
        for (int k = 0; k < v.getDepth(); ++k) for (int j = 0; j < v.getHeight(); ++j) for (int i = 0; i < v.getWidth(); ++i) {
            if (v_weight(i, j, k) > 1e-9) v(i, j, k) /= v_weight(i, j, k);
        }
        for (int k = 0; k < w.getDepth(); ++k) for (int j = 0; j < w.getHeight(); ++j) for (int i = 0; i < w.getWidth(); ++i) {
            if (w_weight(i, j, k) > 1e-9) w(i, j, k) /= w_weight(i, j, k);
        }
    }
    void GridToParticle(MACGrid& grid, const Grid<float>& u_old, const Grid<float>& v_old, const Grid<float>& w_old, float alpha){
        const float dx = grid.getDx();
        // 1. 预先计算速度变化量网格 (只计算一次，不放在粒子循环里)
        Grid<float> delta_u = grid.u() - u_old;
        Grid<float> delta_v = grid.v() - v_old;
        Grid<float> delta_w = grid.w() - w_old;

        // 2. 遍历粒子应用混合公式
        for (auto& p : grid.particles()){
            // PIC 部分：直接从当前网格插值
            Vector3f vel_PIC = Advector::get_velocity_at(grid, p.position);

            // FLIP 部分：使用新重载函数，直接传入 delta 网格，没有任何内存分配
            Vector3f delta_vel_interp = Advector::get_velocity_at(delta_u, delta_v, delta_w, dx, p.position);
            Vector3f vel_FLIP = p.velocity + delta_vel_interp;

            // 应用混合公式
            p.velocity = vel_PIC * alpha + vel_FLIP * (1.0f - alpha);
        }
    }
}