# TABI-PB
This the repository of the Fortran code for treecode-accelerated Poisson-Boltzmann Solver 

This software package employs a well-conditioned boundary integral formulation for the electrostatic potential and its normal derivative on the molecular surface. The surface is triangulated and the integral equations are discretized by centroid collocation. The linear system is solved by GMRES iteration and the matrix-vector product is carried out by a Cartesian terraced which reduces the cost from O(N^2) to O(N*logN), where N is the number of faces in the triangulation. The TABI solver can be applied to compute the electrostatic potential on molecular surface and solvation energy. 

REFERENCE: 
W. Geng and R. Krasny, A treecode-accelerated boundary integral Poisson-Boltzmann solver for continuum electrostatics of solvated biomolecules, J. Comput. Phys., 247, 62-87 (2013).

J. Chen and W. Geng, On preconditioning the treecode-accelerated boundary integral (TABI) Poisson-Boltzmann solver, J. Comput. Phys., 373, 750-762 (2018).
