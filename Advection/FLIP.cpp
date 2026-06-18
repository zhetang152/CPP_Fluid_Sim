#include "FLIP.h"

#include <cmath>

#include "advection.h"

namespace {
    void splatComponent(const Point3f& grid_pos, float value, Grid<float>& accum, Grid<float>& weights) {
        const int base_i = static_cast<int>(std::floor(grid_pos.x));
        const int base_j = static_cast<int>(std::floor(grid_pos.y));
        const int base_k = static_cast<int>(std::floor(grid_pos.z));

        const float fx = grid_pos.x - base_i;
        const float fy = grid_pos.y - base_j;
        const float fz = grid_pos.z - base_k;

        for (int dk = 0; dk <= 1; ++dk) {
            for (int dj = 0; dj <= 1; ++dj) {
                for (int di = 0; di <= 1; ++di) {
                    const int i = base_i + di;
                    const int j = base_j + dj;
                    const int k = base_k + dk;
                    if (i < 0 || i >= accum.getWidth() ||
                        j < 0 || j >= accum.getHeight() ||
                        k < 0 || k >= accum.getDepth()) {
                        continue;
                    }

                    const float wx = (di == 0) ? 1.0f - fx : fx;
                    const float wy = (dj == 0) ? 1.0f - fy : fy;
                    const float wz = (dk == 0) ? 1.0f - fz : fz;
                    const float weight = wx * wy * wz;

                    accum(i, j, k) += value * weight;
                    weights(i, j, k) += weight;
                }
            }
        }
    }

    void normalizeComponent(Grid<float>& values, const Grid<float>& weights) {
        for (int k = 0; k < values.getDepth(); ++k) {
            for (int j = 0; j < values.getHeight(); ++j) {
                for (int i = 0; i < values.getWidth(); ++i) {
                    if (weights(i, j, k) > 1e-9f) {
                        values(i, j, k) /= weights(i, j, k);
                    }
                }
            }
        }
    }
}

namespace FLIPSolver {
    void ParticleToGrid(MACGrid& grid) {
        auto& u = grid.u();
        auto& v = grid.v();
        auto& w = grid.w();

        u.fill(0.0f);
        v.fill(0.0f);
        w.fill(0.0f);

        Grid<float> u_weight(u.getWidth(), u.getHeight(), u.getDepth(), 0.0f);
        Grid<float> v_weight(v.getWidth(), v.getHeight(), v.getDepth(), 0.0f);
        Grid<float> w_weight(w.getWidth(), w.getHeight(), w.getDepth(), 0.0f);

        const float dx = grid.getDx();
        for (const auto& particle : grid.particles()) {
            const Point3f p = particle.position / dx;
            splatComponent(Point3f(p.x, p.y - 0.5f, p.z - 0.5f), particle.velocity.x, u, u_weight);
            splatComponent(Point3f(p.x - 0.5f, p.y, p.z - 0.5f), particle.velocity.y, v, v_weight);
            splatComponent(Point3f(p.x - 0.5f, p.y - 0.5f, p.z), particle.velocity.z, w, w_weight);
        }

        normalizeComponent(u, u_weight);
        normalizeComponent(v, v_weight);
        normalizeComponent(w, w_weight);
    }

    void GridToParticle(
        MACGrid& grid,
        const Grid<float>& u_old,
        const Grid<float>& v_old,
        const Grid<float>& w_old,
        float alpha) {
        const float dx = grid.getDx();
        Grid<float> delta_u = grid.u() - u_old;
        Grid<float> delta_v = grid.v() - v_old;
        Grid<float> delta_w = grid.w() - w_old;

        for (auto& p : grid.particles()) {
            Vector3f vel_PIC = Advector::get_velocity_at(grid, p.position);
            Vector3f delta_vel_interp = Advector::get_velocity_at(delta_u, delta_v, delta_w, dx, p.position);
            Vector3f vel_FLIP = p.velocity + delta_vel_interp;
            p.velocity = vel_PIC * alpha + vel_FLIP * (1.0f - alpha);
        }
    }
}
