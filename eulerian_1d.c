#include "shared_euler.h"
#include "eulerian_1d.h"
#include "array_handling.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "io.h"


void init_sod_grid_1d(ConservedStateVector *U_CELL,
                        ConservedStateVector UL , ConservedStateVector UR ,
                        int gridSize,double x_min,double x_max)
{
    //generalized dx now based on range of x
    const double dx = (x_max - x_min) / gridSize;

    for (int i = 0; i < gridSize; ++i)
    {
        //center for cell xi - now based on range of x
        const double xi = x_min + (i + 0.5 ) * dx;

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


void init_shu_osher_1d(ConservedStateVector *U_CELL,
    ConservedStateVector UL , ConservedStateVector UR,int gridSize,double x_min,double x_max)
{
    //generalized dx now based on range of x
    const double dx = (x_max - x_min) / gridSize;

    for (int i = 0; i < gridSize; ++i)
    {
        //center for cell xi - now based on range of x
        const double xi = x_min + (i + 0.5 ) * dx;

        if (xi < -4)
        {
            U_CELL[i] = UL;
        }
        else
        {

            UR.rho = 1+0.2*sin(5 * xi);
            U_CELL[i] = UR;

        }
    }
}

//part of SSP-RK2 - just the LUi part
void compute_rhs(ConservedStateVector *U_CELL,ConservedStateVector *RHS, PhysicalFluxVector *F_HAT_CELL,int gridSize,double gamma,double dx)
{


    compute_euler_fluxes(U_CELL, F_HAT_CELL, gridSize, gamma);
    for (int i = 1; i <= gridSize - 2; ++i)
    {
        RHS[i].rho =
            -(F_HAT_CELL[i].mass_flux
            - F_HAT_CELL[i-1].mass_flux) / dx;

        RHS[i].momentum =
            -(F_HAT_CELL[i].momentum_flux
            - F_HAT_CELL[i-1].momentum_flux) / dx;

        RHS[i].energy =
            -(F_HAT_CELL[i].energy_flux
            - F_HAT_CELL[i-1].energy_flux) / dx;
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
    int gridSize, double gamma, double t_final, double x_min,double x_max)
{


    const double dx = (x_max - x_min) / gridSize;
    double t = 0;


    ConservedStateVector *U_CELL_NEW = malloc( (gridSize) * sizeof(ConservedStateVector));
    ConservedStateVector *U_STAGE = malloc( (gridSize) * sizeof(ConservedStateVector));
    ConservedStateVector *RHS = malloc( (gridSize) * sizeof(ConservedStateVector));

    // write_state_to_csv(U_CELL,gridSize,gamma,t,dx);
    while (t < t_final)
    {
        double dt = compute_dt_CFL(U_CELL,dx,gridSize,gamma);



        if (t + dt > t_final)
        {
            dt = t_final -t;
        }

        U_STAGE = copy_1d_array(U_CELL,U_STAGE,gridSize);

        //Calcing L(U^n) here
        compute_rhs(U_CELL,RHS,F_HAT_CELL,gridSize,gamma,dx);

        //temp boundary treatment
        //this now uses RHS & new SSP formulation - this is first stage building
        //U^1 here
        for (int i = 1; i <= gridSize - 2; ++i)
        {

            U_STAGE[i].rho = U_CELL[i].rho + dt * RHS[i].rho;

            U_STAGE[i].momentum = U_CELL[i].momentum + dt * RHS[i].momentum;

            U_STAGE[i].energy = U_CELL[i].energy + dt * RHS[i].energy;

        }

        U_CELL_NEW = copy_1d_array(U_CELL,U_CELL_NEW,gridSize);
        //RHS = L(U^(1)
        compute_rhs(U_STAGE,RHS,F_HAT_CELL,gridSize,gamma,dx);

        //building U^(n+1)
        for (int i = 1; i <= gridSize - 2; ++i)
        {
            U_CELL_NEW[i].rho = 0.5 * U_CELL[i].rho + 0.5 * (U_STAGE[i].rho + dt * RHS[i].rho);

            U_CELL_NEW[i].momentum = 0.5 * U_CELL[i].momentum + 0.5 * (U_STAGE[i].momentum + dt * RHS[i].momentum);

            U_CELL_NEW[i].energy = 0.5 * U_CELL[i].energy + 0.5 * (U_STAGE[i].energy + dt * RHS[i].energy);
        }

        U_CELL = copy_1d_array(U_CELL_NEW,U_CELL,gridSize);
        t += dt;

    }
    write_state_to_csv(U_CELL,gridSize,gamma,t,dx,x_min);
    const char *output = "outputs/output.csv";
    printf("Writing to file: %s\n", output);
    plot1d_csv_shu(output,gridSize,t_final);
    // plot1d_csv_sod(output);


    free(U_CELL_NEW);
    free(U_STAGE);
    free(RHS);
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

    //C_cfl = 0.5
    double dt = 0.5 * (dx / a_max);

    return dt;
}