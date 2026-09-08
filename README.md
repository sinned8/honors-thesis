#       Honors Thesis
## A comparative computational study of Lagrangian and Eulerian numerical methods for compressible flow simulation

This project will investigate Lagrangian and Eulerian numerical methods for hypersonic compressible flow simulation as a problem in scientific computing and computational modeling. 
The thesis will involve implementing representative versions of the two approaches in a compiled, typed language and compare them on a one-dimensional modified Shu–Osher shock-density interaction problem at Mach 6. The comparison will focus on algorithmic structure, shock resolution, adaptivity, and computational efficiency, with the goal of analyzing how different computational representations and numerical strategies affect accuracy, flexibility, and performance in simulation. If time permits, the project may also explore differences in computational performance on multicore CPU or GPU-based hardware.

## Methodological background

The thesis is grounded on some of the numerical works from the following papers:

- Gao, T., Liang, T., & Fu, L. (2023). A new smoothed particle hydrodynamics method based on high-order moving-least-square targeted essentially non-oscillatory scheme for compressible flows (2306.00514V1). arXiv. https://doi.org/10.1016/j.jcp.2023.112270
- Natarajan, M., Grout, R., Zhang, W., & Day, M. (2022). A moving embedded boundary approach for the compressible navier-stokes equations in a block-structured adaptive refinement framework. Journal of Computational Physics, 465, 111315. https://doi.org/10.1016/j.jcp.2022.111315


## External Reference Code

The `external/` directory contains third-party code used for verification and comparison purposes.

### Riemann Book

The `riemann_book/` repository was cloned into `external/` and is used to generate the exact Riemann solution for the Sod shock tube. This provides a reference solution for validating the numerical results produced by the Eulerian Rusanov solver.

The external repository is not part of the numerical solver developed for this thesis and is used only for verification and plotting comparisons.

A small number of compatibility changes were required to run the repository with the current Python, NumPy, SciPy, and IPython environment. These changes are documented in the [external compatibility notes](external/README.md).


Source: [Riemann Book repository](<https://github.com/clawpack/riemann_book>)