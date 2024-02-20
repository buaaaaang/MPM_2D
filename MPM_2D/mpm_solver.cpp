//
//  mpm_solver.cpp
//  MPM_2D
//
//  Created by LeeSangheon on 2024/02/14.
//

#include "mpm_solver.hpp"

void Solver::add_object(Vec center, int num, float4 color) {
    for (int i=0; i<num; i++) {
        particles.push_back(Particle(
            (Vec::randomVector()*2.0-Vec(1))*0.08 + center));
    }
    colors.push_back(ParticleColor(color, numParticles, numParticles+num));
    numParticles += num;
}

float4 Solver::getParticleColor(int i) {
    while ( i >= colors[colorFinder].end ) colorFinder++;
    while ( i < colors[colorFinder].start ) colorFinder--;
    return colors[colorFinder].color;
}

void Solver::advance() {
    std::memset(grid, 0, sizeof(grid));
    
    for (auto &p : particles) {
        VectorND<2, int> base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();
        Vec fx = p.x * inv_dx - base_coord.cast<float>();
        Vec w[3] = {
            0.5f * (Vec(1.5) - fx).sqr(),
            Vec(0.75) - (fx - Vec(1.0)).sqr(),
            0.5f * (fx - Vec(0.5)).sqr()
        };
        auto e = std::exp(hardening * (1.0f - p.Jp));
        auto mu = mu_0 * e;
        auto lambda = lambda_0 * e;
        float J = determinant(p.F);
        Mat r, s;
        polar_decomp(p.F, r, s);
        float Dinv = 4 * inv_dx * inv_dx;
        auto PF = (2 * mu * (p.F-r) * p.F.transposed() + lambda * (J-1) * J);
        auto stress = - (dt * vol) * (Dinv * PF);
        auto affine = stress + particle_mass * p.C;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto dpos = (Vec(i, j) - fx) * Vec(dx);
                VectorND<3, float> mass_x_velocity(p.v * particle_mass, particle_mass);
                grid[base_coord.x + i][base_coord.y + j] +=
                    (w[i].x*w[j].y * (mass_x_velocity + VectorND<3, float>(affine * dpos, 0)));
            }
        }
    }
    
    for(int i = 0; i <= n; i++) {
        for(int j = 0; j <= n; j++) {
            auto &g = grid[i][j];
            if (g[2] > 0) {
                g /= g[2];
                g += dt * VectorND<3, float>(0, -200, 0);
                float boundary = 0.05;
                float x = (float) i / n;
                float y = float(j) / n;
                if (x < boundary || x > 1-boundary || y > 1-boundary) { g = VectorND<3, float>(0); }
                if (y < boundary) { g[1] = std::max(0.0f, g[1]); }
            }
        }
    }

    for (auto &p : particles) {
        VectorND<2, int> base_coord = (p.x * inv_dx - Vec(0.5f)).cast<int>();
        Vec fx = p.x * inv_dx - base_coord.cast<float>();
        Vec w[3] = { Vec(0.5) * (Vec(1.5) - fx).sqr(),
                     Vec(0.75) - (fx - Vec(1.0)).sqr(),
                     Vec(0.5) * (fx - Vec(0.5)).sqr() };
        p.C = Mat(0);
        p.v = Vec(0);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                auto dpos = (Vec(i, j) - fx);
                auto grid_v = Vec(grid[base_coord.x + i][base_coord.y + j]);
                auto weight = w[i].x * w[j].y;
                p.v += weight * grid_v;
                p.C += 4 * inv_dx * Mat::outer_product(weight * grid_v, dpos);
            }
        }
        p.x += dt * p.v;
        auto F = (Mat(1) + dt * p.C) * p.F;
        Mat svd_u, sig, svd_v;
        svd(F, svd_u, sig, svd_v);
        for (int i = 0; i < 2 * int(plastic); i++) {
          sig[i][i] = clamp(sig[i][i], 1.0f - 2.5e-2f, 1.0f + 7.5e-3f);
        }
        float oldJ = determinant(F);
        F = svd_u * sig * svd_v.transposed();
        float Jp_new = clamp(p.Jp * oldJ / determinant(F), 0.6f, 20.0f);
        p.Jp = Jp_new;
        p.F = F;
      }
}

void Solver::nextFrame() {
    while (frame_time+frame_dt>time) {
        advance();
        time += dt;
    }
    frame_time += frame_dt;
}
