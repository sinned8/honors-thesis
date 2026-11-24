#include "shared_euler.h"
#include "eulerian_1d.h"
#include "lagrangian_1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{
    double gamma = 1.4;
    int gridSize = 1600;
    //introduce ghost cells for WENO3 since out original indexing had 0 and gsize-1 as physical cells
    // that remained unchanged... instead cells 0,1 & N+2,N+3 are L and R ghost cells
    int totalSize = gridSize + 2 * NGHOST;
    double x_min_shu = -5;
    double x_max_shu = 5;
    double x_min_sod = 0;
    double x_max_sod = 1;


    //***********************************************************|
    //              SOD TEST & SOD ICs                           |
    //***********************************************************|
    //
    // PrimitiveStateVector W1 = {1.0, 0, 1.0};
    // PrimitiveStateVector W2 = {0.125, 0, 0.1};
    //
    // ConservedStateVector U = prim_to_conserved(W1,gamma);
    // ConservedStateVector U2 = prim_to_conserved(W2,gamma);
    // ConservedStateVector *U_CELL = malloc(totalSize * sizeof(ConservedStateVector));
    // PhysicalFluxVector *F_HAT_CELL = malloc( (totalSize - 1) * sizeof(PhysicalFluxVector));
    //
    // init_sod_grid_1d(U_CELL,U,U2,gridSize,x_min_sod,x_max_sod);
    // compute_euler_fluxes(U_CELL,F_HAT_CELL,gridSize,gamma);
    // update_euler_grid(U_CELL,F_HAT_CELL,gridSize,gamma,0.2,x_min_sod,x_max_sod,false);
    // free(U_CELL);
    // free(F_HAT_CELL);


    //**********************************************************|
    //          SHU OSHER TEST & ICs                            |
    //**********************************************************|

    // PrimitiveStateVector W1 = {3.857143, 2.629369, 10.33333};
    // //rho is converted to init_shu_osher_1d
    // PrimitiveStateVector W2 = {1, 0, 1};
    //
    // ConservedStateVector U = prim_to_conserved(W1,gamma);
    // ConservedStateVector U2 = prim_to_conserved(W2,gamma);
    // ConservedStateVector *U_CELL = malloc(totalSize * sizeof(ConservedStateVector));
    // PhysicalFluxVector *F_HAT_CELL = malloc( (totalSize - 1) * sizeof(PhysicalFluxVector));
    //
    // init_shu_osher_1d(U_CELL,U,U2,gridSize,x_min_shu,x_max_shu);
    // update_euler_grid(U_CELL,F_HAT_CELL,gridSize,gamma,1.8,x_min_shu,x_max_shu,false);
    // free(U_CELL);
    // free(F_HAT_CELL);

    // double out1 = weno3_left(1,1,10);
    // double out2 = weno3_right(10,1,1);
    // printf("%lf",out2);
    //
    // // CONSTANT TEST CASE
    // PrimitiveStateVector W1 = {1.0, 1.0, 1.0};
    // PrimitiveStateVector W2 = {1.0, 1.0, 1.0};
    //
    // ConservedStateVector U = prim_to_conserved(W1,gamma);
    // ConservedStateVector U2 = prim_to_conserved(W2,gamma);
    // ConservedStateVector *U_CELL = malloc(totalSize * sizeof(ConservedStateVector));
    // PhysicalFluxVector *F_HAT_CELL = malloc( (totalSize - 1) * sizeof(PhysicalFluxVector));
    //
    // init_sod_grid_1d(U_CELL,U,U2,gridSize,x_min_sod,x_max_sod);
    // update_euler_grid(U_CELL,F_HAT_CELL,gridSize,gamma,0.2,x_min_sod,x_max_sod);
    // free(U_CELL);
    // free(F_HAT_CELL);

    //**********************************************************************|
    //          Kadioglu,Sussman et al.                                     |
    //          Smooth test                                                 |
    // https://www.sciencedirect.com/science/article/pii/S0021999105001798  |
    //**********************************************************************|
    // double x_min_smooth = 0.0;
    // double x_max_smooth = 1.0;
    //
    // int smoothGridSize = 3200;
    // int smoothTotalSize = smoothGridSize + 2 * NGHOST;
    //
    // ConservedStateVector *U_CELL =
    //     malloc(smoothTotalSize * sizeof(ConservedStateVector));
    //
    // PhysicalFluxVector *F_HAT_CELL =
    //     malloc((smoothTotalSize - 1) * sizeof(PhysicalFluxVector));
    //
    // init_smoothtest_1d(U_CELL, smoothGridSize, x_min_smooth,
    //                    x_max_smooth, gamma);
    //
    // //initial sound speed
    // double c0 = sqrt(gamma * 1.0e6 / 1.0e-3);
    //
    // //one acoustic crossing of the domain
    // double t_final = 1.0 / c0;
    //
    // update_euler_grid(U_CELL, F_HAT_CELL, smoothGridSize,
    //     gamma, t_final, x_min_smooth, x_max_smooth, 1);
    //
    // free(U_CELL);
    // free(F_HAT_CELL);
    //


    //**********************************************************|
    //          Particle testing                                |
    //**********************************************************|
    int particleCount = 16;
    ParticleStruct particles[16];

    initial_particle_system(particles,particleCount,0,1);


    // for (int i = 0; i < particleCount; i++)
    // {
    //     printf("Particle %d: x = %lf, V = %lf\n",
    //            i, particles[i].pos, particles[i].volume);
    // }
    // double xi = 0.0, xj = 0.8;
    //
    // double h = 0.4;
    // double W = wendland_c4_1d(xi,xj,h);
    // double derivW = wendland_c4_deriv_1d(xi,xj,h);
    // printf("\n W = %lf \n W deriv = %lf",W,derivW);

    for (int i = 0; i < particleCount; i++)
    {
        // same velo
        // PrimitiveStateVector W = {1.0,1.0,1.0};

        // expansion test
        // PrimitiveStateVector W = {1.0,particles[i].pos,1.0};

        // compression test
        PrimitiveStateVector W = {1.0,-particles[i].pos,1.0};
        particles[i].U = prim_to_conserved(W,gamma);
    }

    double dVdt = compute_volume_RHS(
        particles,8,particleCount);

    printf("dV/dt = %.10lf\n",dVdt);

}
