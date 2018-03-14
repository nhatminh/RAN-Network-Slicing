# RAN Network Slicing

#### This module is built on Julia with JuliaOpt packages (JUMP) for solving optimization.

1. Setting2,3,4,5,6.jl consist traffic generation, and arrival rates of 2,3,4,5,6 tenants, respectively.

2.  EdgeSlicing2,3,4,5,6.jl are the main code files for each scenario.
  * REUSED_INPUT = true (using the pre-generated input, e.g. weights2.h5)
  * FLOT_ALL_FIGS = true (only plot graphs without running optimization)
  * ALG_MODE = 1 #1, 2, 3, 4 (algorithm schemes, we are using 1 (**PBCD**) and 3(**JP-ADMM**) )
  
3.  EdgeSlicing_ExhaustiveSearch.jl is exhaustive search 

4.  Plot_Fig2.jl (for 2 tenants) and Plot_Fig_general.jl (for more than 2 tenants)

  * Folder figs2,3,4,5,6 store figures and output files (e.g. result1.h5 - output from algorithm 1) after running problem.
  
An additional statistical plot via seaborn in python by running `plot.py`

