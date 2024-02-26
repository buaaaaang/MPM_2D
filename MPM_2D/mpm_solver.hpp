//
//  mpm_solver.hpp
//  MPM_2D
//
//  Created by LeeSangheon on 2024/02/14.
//

#ifndef mpm_solver_hpp
#define mpm_solver_hpp

#include <simd/simd.h>

#include <stdio.h>
#include <vector>
#include "utils.hpp"

using Vec = VectorND<2, float>;
using Mat = MatrixND<2, float>;
using float4 = simd::float4;

struct Particle {
  Vec x, v; // Position and velocity
  Mat F; // Deformation gradient
  Mat C; // Affine momentum from APIC
  float Jp; // Determinant of the deformation gradient (i.e. volume)

  Particle(Vec x, Vec v=Vec(0)) :
    x(x), v(v), F(1), C(0), Jp(1) {}
};

struct ParticleColor {
    union {
        float4 color;
        struct {
            float r, g, b, a;
        };
    };
    int start, end;
    
    ParticleColor(float4 color, int start, int end) :
        color(color), start(start), end(end) {}
};

class Solver {
public:
    void advance();
    void nextFrame();
    void add_square(Vec center, int numParticles, float4 color, float size = 0.20);
    void add_circle(Vec center, int numParticles, float4 color, float size = 0.20);
    float getParticleX(int i) { return particles[i].x.x; };
    float getParticleY(int i) { return particles[i].x.y; };
    float4 getParticleColor(int i);
    int numParticles;
    int colorFinder;
    float time;
    float frame_time;
    Solver() : numParticles(0), colorFinder(0), time(0.0), frame_time(0.0) {
        add_square(Vec(0.55,0.3), 800, { 0.906, 0.337, 0.239, 1.0 });
        add_square(Vec(0.60,0.8), 800, { 0.91, 0.694, 0.204, 1.0 });
        add_circle(Vec(0.40,0.55), 800, { 0.031, 0.486, 0.502, 1.0 });
    }
private:
    static const int n = 80; // grid resolution
    std::vector<Particle> particles;
    std::vector<ParticleColor> colors;
    const float dt = 1e-4;
    const float frame_dt = 1e-3;
    const float dx = 1.0 / n;
    const float inv_dx = 1.0 / dx;
    const float particle_mass = 1.0;
    const float vol = 1.0;        // Particle Volume
    const float hardening = 10.0; // Snow hardening factor
    const float E = 1e4;          // Young's Modulus
    const float nu = 0.2;         // Poisson ratio
    const bool plastic = false;
    const float mu_0 = E/(2*(1+nu));
    const float lambda_0 = E*nu/((1+nu)*(1-2*nu));
    VectorND<3, float> grid[n+1][n+1];
};

#endif /* mpm_solver_hpp */
