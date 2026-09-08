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