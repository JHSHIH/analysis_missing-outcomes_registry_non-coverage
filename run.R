# Required source files
sf1 <- "/home/jshih/registry/simu/source/source.R"
sf2 <- "/home/jshih/registry/simu/source/Decode_state.R"
sf3 <- "/home/jshih/registry/simu/source/USRT_VPR_funcs.R"
sf4<- "/home/jshih/registry/simu/source/source.so"

source(sf1) # Needed for calling runSims_genSwarmFiles below

# Folder where the job files will be written to
cdir0  <- "/home/jshih/registry/simu/source/call3/"
odir0  <- "/home/jshih/registry/simu/source/out3/"


STR           <- "sim_Biostatistics"
NREP          <- 1
NBOOT         <- 200
NREPPERJOB    <- 1
NBOOTPERJOB   <- 200


sourceFiles   <- c(sf1, sf2, sf3, sf4)

# seeds for replications and bootstrap samples
set.seed(123)
seeds.rep  <- sample(1:999999999, NREP, replace=FALSE)
seeds.boot <- sample(1:999999999, NBOOT, replace=FALSE)
#########################################################################
# simulation parameters
# All parameters have default values. The default values can be seen by calling check_parms(NULL).
# par1 and par2 are used in:
# impute.method.1.cmpr2(par1,…) 
# impute.method.V.cmpr2(par1,…)
# sample.coxph2.cmpr2(par2,…)
# sample.coxph.V2.cmpr2(par2,…) 

#Moving-stayer model parameters
ms.pars      <- list()
ms.pars[[1]] <- c(0.39,0.97,0.98,0.2,0.4)
#ms.pars[[2]] <- c(0.39,0.9,0.9,0.2,0.4)
#ms.pars[[3]] <- ms.pars[[1]] 
#ms.pars[[4]] <- ms.pars[[2]] 

####Parameter specifications ### ### ###
#Unused parameters: a, a2,b, b2,me.cen,me.coef
#beta: Cox regression coefficients of interest
#epsilon: tolerance level for convergence 
#lambda, lambda2: constant hazard for failure time of type 1 and 2, respectively
#M: # of imputations for missing residence data when par.impute=T
#M2: # of imputations for missing registry-identified failure time 
#MAXT: 300 months (25 years) of followup
#me.sd: level of measurement error
#par.impute: if par.impute=F, impute missing residence status by posterior mode; if par.impute=T, impute missing residence #status from the posterior distribution (M times)
#p_00,p_11,pi_0,pi_1,r_0: initial values of the moving stayer model parameter estimates
#par1: argument passed to impute.1.method.cmpr2
#       if par1=NULL, impute both failure time and failure type if self-reported failure time is censored. Otherwise, impute only failure time and failure type of imputed failure time is the same as the failure type of self-reported failure time.
#       if par1 !=NULL, impute both failure time and failure type.
      
#par2: argument passed to sample.coxph2.cmpr2
#       if par2<0, impute both failure type and failure time 
#       if par2=NULL, impute both failure time and failure type if self-reported failure time is censored. Otherwise impute only failuretime and failure type of imputed failure time is the same as the failure type of self-reported failure time.

#sen, spec: sensitivity and specificity with respect to failure type (1 vs 2) of gold standard failure time and error-prone failure time

### ### #### ### #### #### #### #### #### #### 

parms <- list(a          = 1, 
              a2         = 1/0.22, 
              b          = exp(9), 
              b2         = exp(5),
              beta       = t(matrix(c(0,0,0.2,-0.2,0.4,-0.4),nr=2,nc=3)),
              epsilon    = 0.0000001,
              ib         = 2,
              lambda     = 0.0016, 
              lambda2    = 0.0016*0.75, 
              M          = 1, 
              M2         = 20, 
              MAXT       = 300,
              me.cen     = t(matrix(c(-0.05,-0.075,-0.5,-0.5,-1,-1),nr=2,nc=3)),
              me.coef    = c(0,0.1), 
              me.sd      = c(1), 
              ms.pars    = ms.pars,
              n          = 10000, 
              n.years    = 25,
              p_00       = 0.9712229, 
              p_11       = 0.979527, 
              par1       = NULL,
              par2       = NULL,
              par.impute = F,
              pi_0       = 0.2201757, 
              pi_1       = 0.5087876,
              r_0        = 0.3932528,
              sen        = 0.99, 
              spec       = 0.99)



  op <- list(boot.summary=1)
  runSims_genSwarmFiles(cdir0, odir0, sourceFiles, parms, seeds.rep, seeds.boot=seeds.boot, 
                        nRepPerJob=NREPPERJOB, nBootPerJob=NBOOTPERJOB, op=op, nswarm=1, append=1) 

# To submit jobs, cd to jobfile.dir and enter the command "swarm -f Rjobs1_ -p 2 --time xx:yy:zz -g 2"

# To combine output files, 
# dirs <-list.files(odir0,full.names=TRUE)
# combineSimFiles2(dirs,paste0(odir0,STR),nRepPerFile=1)


