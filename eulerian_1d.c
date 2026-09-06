#include "shared_euler.h"
#include "eulerian_1d.h"
#include "array_handling.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>




void init_euler_grid_1d(ConservedStateVector *U_CELL,
    ConservedStateVector UL , ConservedStateVector UR ,
    int gridSize)
{
    //assuming 0 <= x <= 1
    double dx = 1.0 / gridSize;

    for (int i = 0; i < gridSize; ++i)
    {
        //center for cell xi
        double xi = (i + 0.5) * dx;

        if (xi < 0.5)
        {
            U_CELL[i] = UL;
        }
        else
        {
            U_CELL[i] = UR;
        }
    }

}

void compute_euler_fluxes(ConservedStateVector *U_CELL, PhysicalFluxVector *F_HAT_CELL,
    int gridSize, double gamma)
{

    //computes fluxes between neighboring cells
    for (int i = 0; i < gridSize - 1; ++i)
    {
        F_HAT_CELL[i] = rusanov_flux(U_CELL[i],U_CELL[i+1],gamma);
    }
}

void update_euler_grid(ConservedStateVector *U_CELL, PhysicalFluxVector *F_HAT_CELL,
    int gridSize, double gamma, double t_final)
{

    //assuming 0 <= x <= 1
    const double dx = 1.0 / gridSize;

    double t = 0;

    ConservedStateVector *U_CELL_OLD = malloc( (gridSize) * sizeof(ConservedStateVector));
    ConservedStateVector *U_CELL_NEW = malloc( (gridSize) * sizeof(ConservedStateVector));
    U_CELL_NEW = copy_1d_array(U_CELL,U_CELL_NEW,gridSize);


    while (t < t_final)
    {
         double dt = compute_dt_CFL(U_CELL,dx,gridSize,gamma);
        if (t + dt > t_final)
        {
            dt = t_final -t;
        }
        printf("\ndt: %f\n",dt);
        U_CELL_OLD = copy_1d_array(U_CELL_NEW,U_CELL_OLD,gridSize);

        //temp boundary treatment
        for (int i = 1; i <= gridSize - 2; ++i)
        {
            U_CELL_NEW[i].rho = U_CELL_OLD[i].rho - (dt/dx *(F_HAT_CELL[i].mass_flux - F_HAT_CELL[i-1].mass_flux));
            U_CELL_NEW[i].momentum = U_CELL_OLD[i].momentum - (dt/dx * (F_HAT_CELL[i].momentum_flux - F_HAT_CELL[i-1].momentum_flux));
            U_CELL_NEW[i].energy = U_CELL_OLD[i].energy - (dt/dx * (F_HAT_CELL[i].energy_flux - F_HAT_CELL[i-1].energy_flux));
        }
        U_CELL = copy_1d_array(U_CELL_NEW,U_CELL,gridSize);
        compute_euler_fluxes(U_CELL,F_HAT_CELL,gridSize,gamma);
        t += dt;
    }


    free(U_CELL_OLD);
    free(U_CELL_NEW);
}

double compute_dt_CFL(ConservedStateVector *U_CELL,double dx,int gridSize, double gamma)
{
    double a_max = 0.0;
    for (int i = 0; i < gridSize; ++i)
    {
        PrimitiveStateVector W = conserved_to_prim(U_CELL[i],gamma);
        double c = speed_sound(U_CELL[i],gamma);
        double wave_speed = fabs(W.velocity) + c;

        if (wave_speed > a_max)
        {
            a_max = wave_speed;
        }
    }


    double dt = 0.5 * (dx / a_max);

    return dt;
}