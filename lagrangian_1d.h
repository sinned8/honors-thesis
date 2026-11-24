#ifndef LAGRANGIAN_1D_H
#define LAGRANGIAN_1D_H
#include "shared_euler.h"


void initial_particle_system(ParticleStruct *U,int count,double x_min,double x_max);

double wendland_c4_1d(double xi, double xj, double h);

double wendland_c4_deriv_1d(double xi, double xj, double h);

double compute_volume_RHS(ParticleStruct *P ,int i, double count);

#endif