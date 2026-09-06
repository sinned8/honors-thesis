#include "shared_euler.h"
#include "eulerian_1d.h"
#include <stdio.h>
#include <stdlib.h>


int main()
{
    double gamma = 1.4;
    int gridSize = 100;

    PrimitiveStateVector W1 = {1.0, 0, 1.0};
    PrimitiveStateVector W2 = {0.125, 0, 0.1};

    ConservedStateVector U = prim_to_conserved(W1,gamma);
    ConservedStateVector U2 = prim_to_conserved(W2,gamma);

    // PrimitiveStateVector W2 = conserved_to_prim(U,gamma);
    // printf("rho: %f "
    //        "\nmomentum: %f"
    //        "\nenergy: %f" , U.rho, U.momentum,U.energy );
    // printf("\n\nW2 rho: %f "
    //        "\nW2 velo: %f"
    //        "\nW2 pressure: %f" , W2.rho, W2.velocity,W2.pressure );
    //
    // PhysicalFluxVector F1 = physical_flux(U,gamma);
    //
    // printf("\n\nF1 mass flux: %f "
    // "\nF1 momentum flux : %f"
    // "\nF1 energy flux: %f" , F1.mass_flux, F1.momentum_flux, F1.energy_flux );
    //
    //
    // double c = speed_sound(U,gamma);
    // printf("\n\nSpeed of sound: %f" ,c );


    // PhysicalFluxVector F_HAT = rusanov_flux(U,U2,gamma);
    // printf("\n\nF_HAT mass flux: %f "
    // "\nF_HAT momentum flux : %f"
    // "\nF_HAT energy flux: %f \n" , F_HAT.mass_flux, F_HAT.momentum_flux, F_HAT.energy_flux );

    ConservedStateVector *U_CELL = malloc(gridSize * sizeof(ConservedStateVector));
    PhysicalFluxVector *F_HAT_CELL = malloc( (gridSize - 1) * sizeof(PhysicalFluxVector));

    init_euler_grid_1d(U_CELL,U,U2,gridSize);
    compute_euler_fluxes(U_CELL,F_HAT_CELL,gridSize,gamma);

    // for (int i = 48; i <=51; ++i)
    // {
    //
    //     printf("U_CELL rho: %f"
    //            "\nU_CELL momentum : %f"
    //            "\nU_CELL energy: %f \n" , U_CELL[i].rho,U_CELL[i].momentum,U_CELL[i].energy);
    //
    //     printf("\nF_HAT_CELL mass flux: %f"
    //            "\nF_HAT_CELL momentum flux : %f"
    //            "\nF_HAT_CELL energy flux: %f \n" , F_HAT_CELL[i].mass_flux,F_HAT_CELL[i].momentum_flux,F_HAT_CELL[i].energy_flux);
    // }

    update_euler_grid(U_CELL,F_HAT_CELL,gridSize,gamma,0.2);

    for (int i = 48; i <=51; ++i)
    {

        printf("\nU_CELL rho: %f"
               "\nU_CELL momentum : %f"
               "\nU_CELL energy: %f \n" , U_CELL[i].rho,U_CELL[i].momentum,U_CELL[i].energy);

        printf("\nF_HAT_CELL mass flux: %f"
               "\nF_HAT_CELL momentum flux : %f"
               "\nF_HAT_CELL energy flux: %f \n" , F_HAT_CELL[i].mass_flux,F_HAT_CELL[i].momentum_flux,F_HAT_CELL[i].energy_flux);
    }

    PrimitiveStateVector W3 =
    conserved_to_prim(U_CELL[50], gamma);

    printf("rho = %f\n", W3.rho);
    printf("velocity = %f\n", W3.velocity);
    printf("pressure = %f\n", W3.pressure);

    free(U_CELL);
    free(F_HAT_CELL);
}
