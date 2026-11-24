#include "lagrangian_1d.h"
#include "shared_euler.h"
#include "array_handling.h"
#include "io.h"
#include <math.h>
#include <stdio.h>


void initial_particle_system(ParticleStruct *U,int count,double x_min,double x_max)
{
    double dx = (x_max - x_min) / count;
    for (int i=0;i<count;i++)
    {
        U[i].pos = x_min + (i + 0.5)*dx;
        U[i].volume= dx;
    }
}
double wendland_c4_1d(double xi, double xj, double h)
{
    //normalized distance
    double q  = fabs(xi - xj)/h;
    if (q < 2 )
    {
        double kernel = 3/(4*h) * pow((1-q/2),5) * (2*q*q + 5.0/2.0 * q + 1);
        return kernel;
    }
    else
    {
        return 0;
    }



}

double wendland_c4_deriv_1d(double xi, double xj, double h)
{
    double r = xi - xj;
    double q = fabs(r)/h;

    // particle is outside the kernel support or the kernel
    // deriv is 0 at center also avoid div by 0
    if (q >= 2.0 || r == 0.0)
    {
        return 0.0;
    }
    double distance = r / fabs(r);

    double derivative = -21.0/(8.0*h*h)*(q*(2.0*q+1.0)*pow((1.0-q/2.0),4)* distance);
    return derivative;
}

double compute_volume_RHS(ParticleStruct *P ,int i, double count)
{
    double ui = P[i].U.momentum / P[i].U.rho;
    double Vi = P[i].volume;

    double h = 2 *Vi;

    double dVdt = 0.0;
    for (int j=0;j<count;j++)
    {
        if (i == j)
        {
            continue;
        }
        double distance = fabs(P[i].pos - P[j].pos);
        if (distance < 2*h)
        {
            double uj = P[j].U.momentum / P[j].U.rho;
            double Vj = P[j].volume;
            double dW = wendland_c4_deriv_1d(P[i].pos,P[j].pos,h);

            dVdt += (ui - uj) * dW * Vi * Vj;
        }
    }
    return dVdt;
}
