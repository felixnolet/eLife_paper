# eLife_paper
Code files for eLife paper

General remark: Fortran codes are used to obtain raw data of PDE simulations. Data is analyzed and visualized in Mathematica.

nuc1d.f95:
  Code for nuclear import model in 1d. Used for Fig 3. Positions of nuclei can be changed in the code.

nuc2d.f95:
  Code for nuclear import model in 2d (Fig 3).
  
nuc_freqdependent_1d.f95:
  Code for nuclear import model where the frequency depends on the local concentration. Model in 1d, used for Fig 4.

fhn1d_hpc.f95:
  Fortran code for solving PDE system of FHN model. Gives .txt as output. Used for Fig 4 supp 1+2

cell1d_hpc.f95:
  Fortran code for solving PDE system of CCO model. Gives .txt as output. Used for Fig 4 supp 1
  
FHN_ODEsolver.nb:
  Mathematica code to solve ODE sytem, used for Fig 4 supp 1. Produces a figure of time series of the solutions. Same result can be obtained with fhn1d_hpc.f95 when setting D=0.
