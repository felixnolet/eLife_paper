# eLife_paper
Code files for eLife paper

General remark: Fortran codes are used to obtain raw data of PDE simulations. Data is analyzed and visualized in Mathematica.

fhn1d_hpc.f95:
  Fortran code for solving PDE system of FHN model. Gives .txt as output. Used for Fig 2 and 4.

cell1d_hpc.f95:
  Fortran code for solving PDE system of CCO model. Gives .txt as output. Used for Fig 2 and 4.
  
FHN_ODEsolver.nb:
  Mathematica code to solve ODE sytem, used for Fig 1. Produces a figure of time series of the solutions. Same result can be obtained with fhn1d_hpc.f95 when setting D=0.

toy1d1.f95:
  Code for nuclear import model in 1d. Used for Fig 3. Positions of nuclei can be changed in the code.

toy2d.f95:
  Code for nuclear import model in 2d (Fig 3).
