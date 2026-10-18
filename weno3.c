#include "weno3.h"
#include <math.h>

double weno3_left(double qm1, double q0, double qp1)
{
    //Jiang–Shu set eps to 10^-6 but Schmidtmann et al discuss tying eps to depend on the spacial discretization del x
    // will perhaps investigate in the future about eps fine tuning...
    const double eps = 1e-6;

    double p0 = -0.5 * qm1 + 1.5 * q0;
    double p1 = 0.5 * q0 + 0.5 * qp1;

    double b0 = (q0 - qm1) * (q0 - qm1);
    double b1 = (qp1 - q0) * (qp1 - q0);

    double a0 = (1.0 / 3.0) / ((eps + b0) * (eps + b0));
    double a1 = (2.0 / 3.0) / ((eps + b1) * (eps + b1));

    return (a0 * p0 + a1 * p1) / (a0 + a1);
}

double weno3_right(double q0, double qp1, double qp2)
{
    return weno3_left(qp2, qp1, q0);
}
