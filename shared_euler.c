#include "shared_euler.h"


//switching from primitive states to conserved states
ConservedStateVector prim_to_conserved( PrimitiveStateVector w,double gamma)
{
    ConservedStateVector u;
    u.rho = w.rho;
    u.momentum = w.rho * w.velocity;
    u.energy = (w.pressure / (gamma - 1)) + (0.5 * w.rho * ( w.velocity * w.velocity));

    return u;
}
//switching from conserved states to primitive states
PrimitiveStateVector conserved_to_prim( ConservedStateVector u, double gamma)
{
    PrimitiveStateVector w;
    w.rho = u.rho;
    w.velocity = u.momentum / u.rho;
    w.pressure = (gamma  - 1)*(u.energy - (0.5 * u.rho * (w.velocity * w.velocity) ) );

    return w;
}

//creating a physical flux vector from a conserved state vector u
PhysicalFluxVector physical_flux(ConservedStateVector u, double gamma)
{
    PhysicalFluxVector FLUX_VECTOR;
    PrimitiveStateVector W;
    W = conserved_to_prim(u,gamma);

    FLUX_VECTOR.mass_flux = u.momentum;
    FLUX_VECTOR.momentum_flux = u.rho * (W.velocity * W.velocity) + W.pressure;
    FLUX_VECTOR.energy_flux = W.velocity * (u.energy + W.pressure);

    return FLUX_VECTOR;
}