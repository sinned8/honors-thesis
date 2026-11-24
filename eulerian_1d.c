#include "shared_euler.h"
#include "eulerian_1d.h"
#include "array_handling.h"
#include "weno3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "io.h"


void init_sod_grid_1d(ConservedStateVector *U_CELL,
                        ConservedStateVector UL , ConservedStateVector UR ,
                        int gridSize,double x_min,double x_max)
{
    //generalized dx now based on range of x
    //added proper indexing w/ ghost cells
    const double dx = (x_max - x_min) / gridSize;
    const int totalSize = gridSize + 2 * NGHOST;

    for (int i = 0; i < totalSize; ++i)
    {
        //center for cell xi - now based on range of x and ghost cells
        const double xi = x_min + (i - NGHOST + 0.5) * dx;

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
    //added proper indexing w/ ghost cells
    const double dx = (x_max - x_min) / gridSize;
    const int totalSize = gridSize + 2 * NGHOST;

    for (int i = 0; i < totalSize; ++i)
    {
        //center for cell xi - now based on range of x and ghost cells
        const double xi = x_min + (i - NGHOST + 0.5) * dx;

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

//initializing the smooth test case from Kadioglu,Sussman et al. section 4.1.8
// https://www.sciencedirect.com/science/article/pii/S0021999105001798
void init_smoothtest_1d (ConservedStateVector *U_CELL,int gridSize,double x_min,
    double x_max,double gamma)
{
    const double dx = (x_max - x_min) / gridSize;
    const int totalSize = gridSize + 2 * NGHOST;


    const double pi = acos(-1.0);
    const double p0 = 1.0e6;
    const double rho0 = 1.0e-3;
    const double epsilon = 1.0;

    for (int i = NGHOST; i < gridSize + NGHOST; ++i)
    {
        //center of each physical cell
        const double xi = x_min + (i - NGHOST + 0.5) * dx;

        //smooth pressure perturbation
        const double p1 = 60.0*cos(2.0*pi*xi)
                        + 100.0*sin(4.0*pi*xi);

        const double p = p0 + epsilon*p1;

        //isentropic density relation
        const double rho = rho0 * pow(p / p0, 1.0/gamma);

        PrimitiveStateVector W = {rho, 0.0, p};
        U_CELL[i] = prim_to_conserved(W, gamma);
    }
}

//periodic BC for smooth test
void apply_periodic_ghosts(ConservedStateVector *U_CELL,
                                 int gridSize)
{
    for (int i = 0; i < NGHOST; ++i)
    {
        //left ghosts copy the last physical cells
        U_CELL[i] = U_CELL[gridSize + i];

        //right ghosts copy the first physical cells
        U_CELL[gridSize + NGHOST + i] =
            U_CELL[NGHOST + i];
    }
}



//part of SSP-RK2 - just the LUi part
void compute_rhs(ConservedStateVector *U_CELL,ConservedStateVector *RHS, PhysicalFluxVector *F_HAT_CELL,int gridSize,double gamma,double dx)
{


    compute_euler_fluxes(U_CELL, F_HAT_CELL, gridSize, gamma);
    //now incorporates the ghost cells for WENO3 recon.
    for (int i = NGHOST; i < gridSize + NGHOST; ++i)
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

//checking for valid states - uses ideal-gas pressure relationship for pressure
//returns zero if either the density or pressure < 0 returns 1 otherwise
static int valid_state(ConservedStateVector U, double gamma)
{
    if (!isfinite(U.rho) || !isfinite(U.momentum) ||
        !isfinite(U.energy) || U.rho <= 0.0)
    {
        return 0;
    }

    double pressure = (gamma - 1.0) *
        (U.energy - 0.5 * U.momentum * U.momentum / U.rho);

    if (!isfinite(pressure) || pressure <= 0.0)
    {
        return 0;
    }

    return 1;
}



void compute_euler_fluxes(ConservedStateVector *U_CELL, PhysicalFluxVector *F_HAT_CELL,
    int gridSize, double gamma)
{
   


    //computes fluxes between neighboring cells
    for (int i = NGHOST - 1; i < gridSize + NGHOST; ++i)
    {
        ConservedStateVector UL;
        ConservedStateVector UR;
        //component wise
        UL.rho = weno3_left(U_CELL[i-1].rho,U_CELL[i].rho,U_CELL[i+1].rho);
        UR.rho = weno3_right(U_CELL[i].rho,U_CELL[i+1].rho,U_CELL[i+2].rho);

        UL.momentum = weno3_left(U_CELL[i-1].momentum,U_CELL[i].momentum,U_CELL[i+1].momentum);
        UR.momentum = weno3_right(U_CELL[i].momentum,U_CELL[i+1].momentum,U_CELL[i+2].momentum);

        UL.energy = weno3_left(U_CELL[i-1].energy,U_CELL[i].energy,U_CELL[i+1].energy);
        UR.energy = weno3_right(U_CELL[i].energy,U_CELL[i+1].energy,U_CELL[i+2].energy);

        //checking for valid states here - if density or pressure is < 0
        // our speed of sound func will incur some problems :(
        if (!valid_state(UL, gamma))
        {
            UL = U_CELL[i];
        }
        if (!valid_state(UR, gamma))
        {
            UR = U_CELL[i + 1];
        }


        F_HAT_CELL[i] = rusanov_flux(UL,UR,gamma);
    }
}

void update_euler_grid(ConservedStateVector *U_CELL, PhysicalFluxVector *F_HAT_CELL,
    int gridSize, double gamma, double t_final, double x_min,double x_max,bool periodic)
{


    const double dx = (x_max - x_min) / gridSize;
    const int totalSize = gridSize + 2 * NGHOST;
    double t = 0;


    ConservedStateVector *U_CELL_NEW = malloc( (totalSize) * sizeof(ConservedStateVector));
    ConservedStateVector *U_STAGE = malloc( (totalSize) * sizeof(ConservedStateVector));
    ConservedStateVector *RHS = malloc( (totalSize) * sizeof(ConservedStateVector));

    // write_state_to_csv(U_CELL,gridSize,gamma,t,dx);
    int stepCount  = 0;
    while (t < t_final)
    {
        if (periodic)
        {
            apply_periodic_ghosts(U_CELL,gridSize);
        }
        //done s.t dt only involves actual computational domain
        double dt = compute_dt_CFL(U_CELL+NGHOST,dx,gridSize,gamma);

        if (t + dt > t_final)
        {
            dt = t_final -t;
        }

        U_STAGE = copy_1d_array(U_CELL,U_STAGE,totalSize);

        //Calcing L(U^n) here
        compute_rhs(U_CELL,RHS,F_HAT_CELL,gridSize,gamma,dx);

        //temp boundary treatment
        //this now uses RHS & new SSP formulation - this is first stage building
        //U^1 here
        for (int i = NGHOST; i < gridSize + NGHOST; ++i)
        {

            U_STAGE[i].rho = U_CELL[i].rho + dt * RHS[i].rho;

            U_STAGE[i].momentum = U_CELL[i].momentum + dt * RHS[i].momentum;

            U_STAGE[i].energy = U_CELL[i].energy + dt * RHS[i].energy;

        }

        U_CELL_NEW = copy_1d_array(U_CELL,U_CELL_NEW,totalSize);
        if (periodic)
        {
            apply_periodic_ghosts(U_STAGE,gridSize);
        }

        //RHS = L(U^(1)
        compute_rhs(U_STAGE,RHS,F_HAT_CELL,gridSize,gamma,dx);

        //building U^(n+1)
        for (int i = NGHOST; i < gridSize + NGHOST; ++i)
        {
            U_CELL_NEW[i].rho = 0.5 * U_CELL[i].rho + 0.5 * (U_STAGE[i].rho + dt * RHS[i].rho);

            U_CELL_NEW[i].momentum = 0.5 * U_CELL[i].momentum + 0.5 * (U_STAGE[i].momentum + dt * RHS[i].momentum);

            U_CELL_NEW[i].energy = 0.5 * U_CELL[i].energy + 0.5 * (U_STAGE[i].energy + dt * RHS[i].energy);
        }

        U_CELL = copy_1d_array(U_CELL_NEW,U_CELL,totalSize);
        t += dt;
        stepCount++;

    }
    write_state_to_csv(U_CELL + NGHOST,gridSize,gamma,t,dx,x_min);
    const char *output = "outputs/output.csv";
    printf("Writing to file: %s\n", output);
    // plot1d_csv_shu(output,gridSize,t_final);
    plot1d_csv_sod(output);
    // printf("N = %d, Timesteps = %d, Final time = %.15e\n",
    //    gridSize, stepCount, t);

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