#include "array_handling.h"

ConservedStateVector *copy_1d_array( ConservedStateVector *arraycp, ConservedStateVector*arraypst, int n)
{

    for (int i = 0; i < n; i++)
    {
        arraypst[i] = arraycp[i];
    }
    return arraypst;
}

PrimitiveStateVector *copy_1d_array_prim(ConservedStateVector *arraycp, PrimitiveStateVector*arraypst, int n, double gamma)
{
    for (int i = 0; i < n; i++)
    {
        arraypst[i] = conserved_to_prim(arraycp[i],gamma);
    }
    return arraypst;
}