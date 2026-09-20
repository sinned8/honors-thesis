#ifndef IO_H
#define IO_H
#include "shared_euler.h"

void write_state_to_csv(ConservedStateVector *U_CELL, int gridSize, double gamma,
    double time, double dx, double x_min);


void plot1d_csv_sod(const char * filename);

void plot1d_csv_shu(const char * filename,int N,double t_final);
#endif
