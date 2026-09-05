#include "array_handling.h"

ConservedStateVector *copy_1d_array( ConservedStateVector *arraycp, ConservedStateVector*arraypst, int n)
{


    for (int i = 0; i < n; i++)
    {
        arraypst[i] = arraycp[i];
    }
    return arraypst;
}