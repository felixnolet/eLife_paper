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
  

Table with experiments: 
Overview cell-free extract experiments. For each experiment brief descriptors are provided: date of the experiment, width of the Teflon tube (µm), sperm dilution (nuclei / µL), time/frame (s), length/pixel (um), whether images were obtained by confocal microscopy, specific compound or condition tested in each tube, followed by further information about the concentration of the condition tested and/or the reference of the reagent. Control: extract supplemented with GFP-NLS (~25 µM); MT: microtubule reporter; DNA: purified genomic DNA; STLC: S-Trityl-L-cysteine; IZ: Importazole. Additionally, the presence of cycles and waves in the cell-free extract. C: cycling extract; NC: Non-cycling extract; O: no wave; B: boundary-dricen wave; I: internal pacemaker-driven wave; IB: coexistence of internally and boundary driven waves.
