#include "shared_euler.h"
#include <stdio.h>



int main()
{
    double gamma = 1.4;

    PrimitiveStateVector W1 = {2.0, 3.0, 1.0};

    ConservedStateVector U = prim_to_conserved(W1,gamma);

    PrimitiveStateVector W2 = conserved_to_prim(U,gamma);
    printf("rho: %f "
           "\nmomentum: %f"
           "\nenergy: %f" , U.rho, U.momentum,U.energy );
    printf("\n\nW2 rho: %f "
           "\nW2 velo: %f"
           "\nW2 pressure: %f" , W2.rho, W2.velocity,W2.pressure );

    PhysicalFluxVector F1 = physical_flux(U,gamma);

    printf("\n\nF1 mass flux: %f "
    "\nF1 momentum flux : %f"
    "\nF1 energy flux: %f" , F1.mass_flux, F1.momentum_flux, F1.energy_flux );


    double c = speed_sound(U,gamma);
    printf("\n\nSpeed of sound: %f" ,c );


    PhysicalFluxVector F_HAT = rusanov_flux(U,U,gamma);
    printf("\n\nF_HAT mass flux: %f "
    "\nF_HAT momentum flux : %f"
    "\nF_HAT energy flux: %f" , F_HAT.mass_flux, F_HAT.momentum_flux, F_HAT.energy_flux );


}
