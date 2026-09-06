#ifndef EULERIAN_1D_H
#define EULERIAN_1D_H


void init_euler_grid_1d(ConservedStateVector *U_CELL,
    ConservedStateVector UL , ConservedStateVector UR,int gridSize);

void compute_euler_fluxes(ConservedStateVector *U_CELL, PhysicalFluxVector *F_HAT_CELL,
    int gridSize, double gamma);

void update_euler_grid(ConservedStateVector *U_CELL, PhysicalFluxVector *F_HAT_CELL,
    int gridSize,double gamma,double t_final);

double compute_dt_CFL(ConservedStateVector *U_CELL,double dx, int gridSize,double gamma);
#endif

