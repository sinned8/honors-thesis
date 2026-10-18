#ifndef WENO3_H
#define WENO3_H

double weno3_left(double qm1, double q0, double qp1);
double weno3_right(double q0, double qp1, double qp2);



#endif
