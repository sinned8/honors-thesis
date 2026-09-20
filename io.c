#include "io.h"
#include "shared_euler.h"
#include "array_handling.h"
#include <stdio.h>
#include <stdlib.h>


// functionality just for final state sod test
void write_state_to_csv(ConservedStateVector *U_CELL, int gridSize, double gamma,
    double time, double dx, double x_min)
{
    PrimitiveStateVector *W_CELL = malloc( (gridSize) * sizeof(PrimitiveStateVector));
    W_CELL = copy_1d_array_prim(U_CELL,W_CELL,gridSize,gamma);


    // just generating new csv file
    FILE *fp = fopen("outputs/output.csv", "w");
    if (fp == NULL)
    {
        printf("Error opening file.\n");
        return;
    }
    fprintf(fp, "x,Density,Velocity,Pressure\n");
    for (int i = 0; i < gridSize; i++)
    {
        const double xi = x_min + (i + 0.5 ) * dx;
        fprintf(fp,"%f,%f,%f,%f\n", xi,U_CELL[i].rho,W_CELL[i].velocity,W_CELL[i].pressure);
    }


    fclose(fp);
    free(W_CELL);
}

void plot1d_csv_sod(const char * filename)
{
    char command[512];

    // command to run plot_1d.py and transfer over necessary variables
    snprintf(command, sizeof(command),
             ".\\.venv\\Scripts\\python.exe "
             "\"python scripts\\plot_1d_sod.py\" "
             "\"%s\"",
             filename);

    printf("Running command:\n%s\n", command);

    int result = system(command);

    if (result != 0)
    {
        printf("Error: Python plotting script failed.\n");
    }
}

void plot1d_csv_shu(const char * filename,int N,double t_final)
{
    char command[512];

    // command to run plot_1d.py and transfer over necessary variables
    snprintf(command, sizeof(command),
             ".\\.venv\\Scripts\\python.exe "
             "\"python scripts\\plot_1d_shu.py\" "
             " \"%s\" \"%d\" \"%lf\" ",
             filename,N,t_final);

    printf("Running command:\n%s\n", command);

    int result = system(command);

    if (result != 0)
    {
        printf("Error: Python plotting script failed.\n");
    }
}