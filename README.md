# eLife_paper
Code files for eLife paper

fhn1d_hpc.f95:
  Fortran code for solving PDE system of FHN model. Gives .txt as output. Used for Fig 2 and 4.

cell1d_hpc.f95:
  Fortran code for solving PDE system of CCO model. Gives .txt as output. Used for Fig 2 and 4.
  
FHN_ODEsolver.nb:
  Mathematica code to solve ODE sytem, used for Fig 1. Same result can be obtained with fhn1d_hpc.f95 when setting D=0.

toy1d1.f95:
  Code for nuclear import model in 1d. Used for Fig 3.

toy2d.f95:
  Code for nuclear import model in 2d (Fig 3).
