#include "shared_euler.h"
#include "eulerian_1d.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{
    double gamma = 1.4;
    int gridSize = 200;
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
    // update_euler_grid(U_CELL,F_HAT_CELL,gridSize,gamma,0.2,x_min_sod,x_max_sod);
    //


    //**********************************************************|
    //          SHU OSHER TEST & ICs                            |
    //**********************************************************|

    PrimitiveStateVector W1 = {3.857143, 2.629369, 10.33333};
    //rho is converted to init_shu_osher_1d
    PrimitiveStateVector W2 = {1, 0, 1};

    ConservedStateVector U = prim_to_conserved(W1,gamma);
    ConservedStateVector U2 = prim_to_conserved(W2,gamma);
    ConservedStateVector *U_CELL = malloc(totalSize * sizeof(ConservedStateVector));
    PhysicalFluxVector *F_HAT_CELL = malloc( (totalSize - 1) * sizeof(PhysicalFluxVector));

    init_shu_osher_1d(U_CELL,U,U2,gridSize,x_min_shu,x_max_shu);
    update_euler_grid(U_CELL,F_HAT_CELL,gridSize,gamma,1.8,x_min_shu,x_max_shu);


    // double out1 = weno3_left(1,1,10);
    // double out2 = weno3_right(10,1,1);
    // printf("%lf",out2);

    // CONSTANT TEST CASE
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
    free(U_CELL);
    free(F_HAT_CELL);
}
