#ifndef SHARED_EULER_H
#define SHARED_EULER_H

//conserved state vector
typedef struct {
    double rho;
    double momentum;
    double energy;
} ConservedStateVector;

//primitive state vector
typedef struct {
    double rho;
    double velocity;
    double pressure;
} PrimitiveStateVector;


// physical flux vector
typedef struct {
    double mass_flux;
    double momentum_flux;
    double energy_flux;
} PhysicalFluxVector;


ConservedStateVector prim_to_conserved( PrimitiveStateVector w,double gamma);
PrimitiveStateVector conserved_to_prim( ConservedStateVector u, double gamma);
PhysicalFluxVector physical_flux(ConservedStateVector u, double gamma);

#endif
