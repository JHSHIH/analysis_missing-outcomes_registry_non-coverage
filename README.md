# Imputation_missing-outcomes_registry_non-coverage
An imputation approach for time-to-event analysis subject to missing outcomes due to non-coverage in disease registries
Authors:  Joanna H. Shih, Paul S. Albert, Jason Fine, Danping Liu

Source code for running simulations:

A set of R source code listed below was used to run simulations on Linux operating system.
Decode_state.R,
USRT_VPR_funcs.R,
source.R,
source.so,
run.R

The source code “source.so” may need be recompiled by using the command “R CMD SHLIB source.c”.

How to run simulations:

1.	Execute R command: source(“run.R”). In the “run.R” code file, number of replicates (NREP)=1 and number of bootstrap samples (NBOOT)=200.
2.	To submit jobs, cd to cdir0 (specified in “run.R”) and enter the command “swarm -f Rjobs1_ -p 2 –time xx:yy::zz -g 2”

Simulation results for one replicate:

1.	Simulation results are saved in the directory odir0. To combine output files, run the following commands:
“dirs <- list.files(odir0, full.names=TRUE)”
“combineSimFiles2(dirs, paste0(odir0, STR), nRepPerFile=1)”

 combineSimFiles2 is a function listed in source.R.

2.	Four output files are created after running combineSimFiles2:
File “sim_Biostatistics_DUMY1” contains beta estimates in Cox regression: beta.GS, beta.IR, beta.Hybrid beta.SR, beta.Step_1, beta.NP, beta.Proposed. 

File “sim_Biostatistics_DUMY2” contains standard errors of beta estimates without accounting for variabilities in the mover-stayer model, imputed coverage status calculated from the posterior model and semi-parametric working imputation model.

File “sim_Biostatistics_DUMY3” contains the estimates of the parameters in the mover-stayer model.

File “sim_Biostatistics_DUMY4_boot” contains the bootstrap standard errors of the parameter estimates in the mover-stayer model.

File “sim_Biostatistics_DUMY5_boot” contains the bootstrap standard errors of the beta estimates in Cox regression.

