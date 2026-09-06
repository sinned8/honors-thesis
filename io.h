#ifndef IO_H
#define IO_H
#include "shared_euler.h"

void write_state_to_csv(ConservedStateVector *U_CELL, int gridSize, double gamma, double time,double dx);


void plot1d_csv(const char * filename);
#endif
