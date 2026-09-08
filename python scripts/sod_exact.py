"""
Exact Sod shock-tube reference solution.

Uses the exact Euler Riemann solver provided by:
Ketcheson, LeVeque, and del Razo,
Riemann Problems and Jupyter Solutions.

https://faculty.washington.edu/rjl/riemann_book/

Used only as a reference solution for verification of the
numerical finite-volume solver.
"""

import sys
from pathlib import Path

import numpy as np


# Clawpack / UW Riemann-book reference implementation
RIEMANN_BOOK = Path(__file__).resolve().parent.parent / "external" / "riemann_book"
sys.path.insert(0, str(RIEMANN_BOOK))

from exact_solvers import euler

def exact_sod_solution(x, t, x0=0.5, gamma=1.4):
    """
    Evaluate the exact Sod shock tube solution using the
    Clawpack Riemann Book exact Euler solver.

    Reference:
    Ketcheson, LeVeque, and del Razo,
    Riemann Problems and Jupyter Solutions.
    """

    # Standard Sod primitive states: (density, velocity, pressure)
    rho_l = 1.0
    u_l = 0.0
    p_l = 1.0

    rho_r = 0.125
    u_r = 0.0
    p_r = 0.1

    # Clawpack's exact solver expects conserved states
    q_l = euler.primitive_to_conservative(
        rho_l, u_l, p_l, gamma=gamma
    )
    q_r = euler.primitive_to_conservative(
        rho_r, u_r, p_r, gamma=gamma
    )

    # Construct exact Riemann solution
    states, speeds, reval, wave_types = euler.exact_riemann_solution(
        q_l, q_r, gamma=gamma
    )

    # Riemann solution is expressed in self-similar coordinate xi.
    xi = (x - x0) / t

    # Evaluate exact conserved state at every x
    q_exact = reval(xi)

    # Convert back to rho, u, p for plotting
    rho_exact, u_exact, p_exact = euler.conservative_to_primitive(
        q_exact[0],
        q_exact[1],
        q_exact[2],
        gamma=gamma
    )

    return rho_exact, u_exact, p_exact