# FAWave1D
A set of MATLAB functions and scripts for performing Fourier analysis of the 1D Wave equation. The available discretization methods are Discontinuous Galerkin (DG) (linear/warped elements), Finite Difference (FD), Compact Difference (CD), and Dispersion Relation Preserving (DRP) schemes and Runge Kutta (RK) for time integration. It also includes explicit Bogey and Bailly filters in addition to Pade' filters. 

This version works with warped non-straight elements but in 1D. Which means each element has extra nodes not just two to define the geometry.
