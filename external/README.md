# External Dependencies

## Riemann Problems and Jupyter Solutions

This project uses the exact Euler Riemann solver from:

Ketcheson, LeVeque, and del Razo,  
*Riemann Problems and Jupyter Solutions*

Repository:
https://github.com/clawpack/riemann_book

Website:
https://faculty.washington.edu/rjl/riemann_book/

The solver is used only as an external reference solution for verification of the numerical Sod shock-tube implementation. It is not part of the Eulerian or Lagrangian numerical methods developed for this thesis.

Clone the repository into this directory with:

```bash
git clone https://github.com/clawpack/riemann_book.git riemann_book
```

## Compatibility Changes

### 1. IPython Display Import

In `riemann_book/utils/nbtools.py`, replace:

```python
from IPython.core.display import display
```

with:

```python
from IPython.display import display
```

This change was needed for compatibility with the current IPython version used in this project environment.

### 2. Scalar Conversion in `euler.py`

In `riemann_book/exact_solvers/euler.py`, add:

```python
p = np.asarray(p).item()
```

immediately before:

```python
u = phi_l(p)
```

This change was needed because `scipy.optimize.fsolve` returns a one-element NumPy array, while the remainder of the solver expects a scalar. Without this conversion, the solver raised a NumPy compatibility error when assigning wave speeds.

## Shu–Osher WENO5 Reference Solution

This project also uses an externally generated, high-resolution numerical solution of the Shu–Osher problem from:

Repository:
https://github.com/sulli72/1D_EULER

Reference file:
`SHU_RHO_PROFILE_EXACT.mat`

The reference was generated using a characteristic-wise WENO5 method with 4,000 grid points and is used to compare density profiles at `t = 1.8`. Although the file is named `EXACT`, it is a **numerical reference solution**, not an analytical exact solution.

To reproduce the comparison, download `SHU_RHO_PROFILE_EXACT.mat` from the repository and place it in this `external/` directory. The local analysis scripts use:

- `python scripts/plot_shu_comparison.py` to overlay the reference density profile with the Eulerian solver output.
- `python scripts/shu_convergence.py` to calculate density error norms and observed convergence orders at different grid resolutions.

The MATLAB file contains reference coordinate and density arrays (`xv` and `r`). Since the reference is supplied at discrete points while the finite-volume solver stores cell averages, the convergence script estimates reference cell averages by integrating a piecewise-linear interpolation of the reference data. These comparisons are approximate and depend on the resolution and accuracy of the reference solution.

The `.mat` file is not committed to this repository; obtain it from its original source using the link above. This external reference is not part of the numerical methods implemented for the thesis.
