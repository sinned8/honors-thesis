#include "shared_euler.h"
#include <math.h>

//switching from primitive states to conserved states
ConservedStateVector prim_to_conserved(const PrimitiveStateVector w, const double gamma)
{
    ConservedStateVector u;
    u.rho = w.rho;
    u.momentum = w.rho * w.velocity;
    u.energy = (w.pressure / (gamma - 1)) + (0.5 * w.rho * ( w.velocity * w.velocity));

    return u;
}
//switching from conserved states to primitive states
PrimitiveStateVector conserved_to_prim(const ConservedStateVector u, const double gamma)
{
    PrimitiveStateVector w;
    w.rho = u.rho;
    w.velocity = u.momentum / u.rho;
    w.pressure = (gamma  - 1)*(u.energy - (0.5 * u.rho * (w.velocity * w.velocity) ) );

    return w;
}

//creating a physical flux vector from a conserved state vector u
PhysicalFluxVector physical_flux(const ConservedStateVector u, const double gamma)
{
    PhysicalFluxVector FLUX_VECTOR;
    const PrimitiveStateVector W = conserved_to_prim(u, gamma);

    FLUX_VECTOR.mass_flux = u.momentum;
    FLUX_VECTOR.momentum_flux = u.rho * (W.velocity * W.velocity) + W.pressure;
    FLUX_VECTOR.energy_flux = W.velocity * (u.energy + W.pressure);

    return FLUX_VECTOR;
}

//speed of sound calc
double speed_sound(const ConservedStateVector u, const double gamma)
{
    const PrimitiveStateVector w = conserved_to_prim(u,gamma);

    const double c = sqrt(gamma * (w.pressure / u.rho));
    return c;
}


//Rusanov formula applied to the conserved UL and UR vectors returns F_HAT which is the Rusanov flux vector
PhysicalFluxVector rusanov_flux(ConservedStateVector UL, ConservedStateVector UR, double gamma)
{
    PhysicalFluxVector FUL = physical_flux(UL,gamma);
    PhysicalFluxVector FUR = physical_flux(UR,gamma);
    PhysicalFluxVector F_HAT;

    double wavespeed_L = fabs(UL.momentum / UL.rho) + speed_sound(UL,gamma);
    double wavespeed_R = fabs(UR.momentum / UR.rho) + speed_sound(UR,gamma);

    double a_max = fmax(wavespeed_L,wavespeed_R);

    //mass
    F_HAT.mass_flux = 0.5 * (FUL.mass_flux + FUR.mass_flux) - 0.5 * (a_max * (UR.rho - UL.rho));
    //momentum
    F_HAT.momentum_flux = 0.5 * (FUL.momentum_flux + FUR.momentum_flux) - 0.5 * (a_max * (UR.momentum - UL.momentum));
    //energy
    F_HAT.energy_flux = 0.5 * (FUL.energy_flux + FUR.energy_flux) - 0.5 * (a_max * (UR.energy - UL.energy));

    return F_HAT;
}