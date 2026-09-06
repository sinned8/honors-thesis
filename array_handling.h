#ifndef ARRAY_HANDLING_H
#define ARRAY_HANDLING_H
#include "shared_euler.h"

ConservedStateVector *copy_1d_array( ConservedStateVector *arraycp, ConservedStateVector*arraypst, int n);

PrimitiveStateVector *copy_1d_array_prim(ConservedStateVector *arraycp, PrimitiveStateVector*arraypst, int n,double gamma);


#endif

