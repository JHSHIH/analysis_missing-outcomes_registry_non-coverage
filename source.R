library(survival)
library(prodlim)
library(splines)
library(plyr)
library(rio)

runSims <- function(out.str, seeds.rep, parms, seeds.boot=NULL, op=NULL) {

  # out.str         Output string that includes the complete path. Ex: /data/wheelerwi/simulation/out/sim1
  # seeds.rep       Numeric vector of replication seeds. Each replication should have its own seed.  
  # parms           List of parameters. See the check_parms function for all parameters and their default values.
  # seeds.boot      Numeric vector of bootstrap seeds. The default is 1:500, ie, 500 bootstraps will be done.

  if (!length(out.str)) stop("ERROR with out.str")

  nrep <- length(seeds.rep)
  if (!length(nrep)) stop("ERROR with seeds.rep")
  nboot <- length(seeds.boot)
  if ((nrep > 1) && nboot) stop("ERROR: not yet implemented")

  op <- default.list(op, c("boot.summary", "append", "out.id.boot"), list(1, 0, 1))

  # Check input parameters and assign default values if needed
  parms  <- check_parms(parms)

  out.rda <- paste0(out.str, "_info.rda")
  info    <- list(seeds.rep=seeds.rep, seeds.boot=seeds.boot, parms=parms, 
                      op=op)
  save(info, file=out.rda)

  # Create output file names and remove them if they exist
  out.files <- paste0(out.str, "_DUMY", 1:5)
  if (!op$append) initOutputFiles(out.files, parms)

  if (nboot) {
    out.files.boot <- paste0(out.str, "_DUMY", 1:5, "_boot")
    if (!op$append) initOutputFiles(out.files.boot, parms)
    out.id.boot <- op$out.id.boot
  }

  t00 <- proc.time()
  for (i in 1:nrep) {
    t0 <- proc.time()

    set.seed(seeds.rep[i])
    cat(paste0("rep = ", i, ", seed = ", seeds.rep[i], "\n"))
    out.id <- paste0("rep", i, "_boot0")

    # Call main function, data will be generated and returned
    data_obj <- try(runSim_main(NULL, out.files, parms, out.id=out.id))
    #data_obj <- runSim_main(NULL, out.files, parms, out.id=out.id)
    if ("try_error" %in% class(data_obj)) next

    # Bootstraps
    if (nboot) {
      for (j in 1:nboot) {
        # For the correct boot numbers in output file
        bootj <- j + out.id.boot - 1 
t00 <- proc.time()
        cat(paste0("boot = ", bootj, ", seed = ", seeds.boot[j], "\n"))
        out.id <- paste0("rep", i, "_boot", bootj)

        # Get bootstrap sample
        rows <- runSim_getBootRows(parms$n, seeds.boot[j])

        # Call main function passing in data returned above
        tmp  <- try(runSim_main(data_obj, out.files.boot, parms, boot.rows=rows, out.id=out.id))
        tmp  <- NULL
        gc()
print(proc.time()-t00)
      }
#       if (op$boot.summary) tmp <- try(getBootSE_file(out.files[1], out.files.boot[1], parms, out.files.boot[5]))
       if (op$boot.summary) {tmp <- try(getBootSE_file(out.files[1], out.files.boot[1], parms, out.files.boot[5]))
              tmp <- try(getBootSE_file2(out.files[3], out.files.boot[3], parms, out.files.boot[4]))
  }

    }

    print(proc.time()-t0)
  }

  print(sessionInfo())
  t11 <- proc.time() - t00
  print(t11)
  NULL
}

initOutputFiles <- function(out.files, parms) {

  tmp    <- file.exists(out.files)
  if (any(tmp)) file.remove(out.files[tmp])
  nc.beta  <- ncol(parms$beta)
  betacols <- paste0("beta.ib.", 1:nc.beta)
  zvec     <- paste0(".", 1:2)
  idcol    <- "rep_boot_ims_ib_ime1_ime2"
  
  out.est <- c(idcol, "ims", "NA.1", 
               betacols, 
               "me.sd.ime1", "me.coef.ime2", 
               paste0("beta.GS", zvec), paste0("beta.IR", zvec), 
               paste0("beta.Hybrid", zvec), paste0("beta.SR", zvec),
               paste0("beta.Step_1", zvec), "NA.2", "NA.3", 
               paste0("beta.NP", zvec), "NA.4", "NA.5", 
               paste0("beta.Proposed", zvec), "NA.6", "NA.7")
  write(out.est, file=out.files[1], sep="\t", ncol=length(out.est), append=T)

  out.se <- c(idcol, "ims", "NA.1",
              betacols,
              "me.sd.ime1", "me.coef.ime2", 
              paste0("se.beta.GS", zvec), paste0("se.beta.IR", zvec),
              paste0("se.beta.Hybrid", zvec), paste0("se.beta.SR", zvec),
              paste0("se.beta.Step_1", zvec), "NA.2", "NA.3", 
              paste0("se.beta.NP", zvec), "NA.4", "NA.5", 
              paste0("se.beta.Proposed", zvec), "NA.6", "NA.7",
              paste0("covp.GS", zvec), paste0("covp.IR", zvec), 
              paste0("covp.Hybrid", zvec), paste0("covp.SR", zvec), 
              paste0("covp.Step_1", zvec), paste0("covp.NP", zvec), 
              paste0("covp.Proposed", zvec))
  write(out.se, file=out.files[2], sep="\t", ncol=length(out.se), append=T)

  #out.par <- c(idcol, "ims", "me.sd.ime1",
  #             "me.coef.ime2", "tmp1", "tmp2",
  #             "tmp3",
  #             "tau1", "tmp4", "tau2", 
  #             "pi0", "pi1", "p0","p1","p00","p11")
  out.par <- c(idcol, "ims", "me.sd.ime1",
                "me.coef.ime2",  "tau1",  "tau2", 
               "pi0", "pi1", "p0","p1","p00","p11")
  write(out.par, file=out.files[3], sep="\t", ncol=length(out.par), append=T)

  NULL
}

runSim_covs.dat <- function(n, beta, ib, lambda, lambda2, year.entry, MAXT) {

  z1=ifelse(runif(n)>0.5,1,-1)
  z2=rnorm(n,0,1)
  #t=-log(runif(n))/(exp(cbind(z1,z2)%*%beta[ib,])*lambda)
  #tc<- -log(runif(n))/lambda2 #time to competing event
  t.all=-log(runif(n))/(exp(cbind(z1,z2)%*%beta[ib,])*lambda+lambda2)
  year.entry<-rep(7,n)
  month.entry<-year.entry*12
  maxt2<-MAXT-month.entry
  #u1=rnorm(n,180,40)
  u1=rnorm(n,150,25)
  u2=maxt2
  c=ifelse(runif(n)<=0.1,u1,u2)
  c<-ifelse(c<maxt2,c,maxt2)
  #x<-ceiling(apply(cbind(t.all,c),1,min))
  x<-apply(cbind(t.all,c),1,min) #integer months are not required
  delta<-ifelse(t.all<c,1,0)
  h.all<-exp(cbind(z1,z2)%*%beta[ib,])*lambda+lambda2
  h.ratio<-exp(cbind(z1,z2)%*%beta[ib,])*lambda/h.all
  type<-ifelse(runif(n)<=h.ratio,1,2)
  delta[delta!=0]<-type[delta!=0]

  covs.dat=data.frame(1:n,rep(NA,n),rep(NA,n),x,delta,z1,z2)
  names(covs.dat)=c("ID","x.obs","delta.obs","x.true","delta.true","z1","z2")

  list(covs.dat=covs.dat, t.all=t.all, z2=z2, 
       year.entry=year.entry, month.entry=month.entry, 
       c=c, type=type)
}

runSim_getBootRows <- function(n, seed) {

  set.seed(seed)
  ret <- sample(1:n, n, replace=TRUE)
  ret

} # END: runSim_getBootRows


runSim_main <- function(data_obj, out.files, parms, boot.rows=NULL, out.id=1) {

  M          <- parms$M
  M2         <- parms$M2
  n          <- parms$n
  beta       <- parms$beta
  lambda     <- parms$lambda
  lambda2    <- parms$lambda2
  a          <- parms$a
  b          <- parms$b #scale
  a2         <- parms$a2
  b2         <- parms$b2
  MAXT       <- parms$MAXT
  me.sd      <- parms$me.sd
  me.coef    <- parms$me.coef
  me.cen     <- parms$me.cen
  ms.pars    <- parms$ms.pars
  n.years    <- parms$n.years
  r_0        <- parms$r_0
  p_00       <- parms$p_00 
  p_11       <- parms$p_11 
  pi_0       <- parms$pi_0
  pi_1       <- parms$pi_1
  epsilon    <- parms$epsilon
  sen        <- parms$sen
  spec       <- parms$spec
  ib         <- parms$ib
  par1       <- parms[["par1", exact=TRUE]]
  par2       <- parms[["par2", exact=TRUE]]
  par.impute <- parms$par.impute
  ims_INDEX  <- parms[["ims_INDEX", exact=TRUE]]
  ime1_INDEX <- parms[["ime1_INDEX", exact=TRUE]]

  if (is.null(data_obj)) {
    bootflag    <- FALSE
  } else {
    bootflag    <- TRUE
  }
  if (!length(ims_INDEX))  ims_INDEX  <- 1:length(ms.pars)
  if (!length(ime1_INDEX)) ime1_INDEX <- 1:length(me.sd)
  if (bootflag) {
    if (length(ims_INDEX) > 1)  stop("ERROR: when bootstrapping, length(ims_INDEX) must be 1")
    if (length(ime1_INDEX) > 1) stop("ERROR: when bootstrapping, length(ime1_INDEX) must be 1")
  }

  irep              <- 1
  out.est           <- NULL
  out.se            <- NULL
  sim_data_obs_list <- list()
  
  # Get the data
  if (!bootflag) {
    # Generate data
    tmp         <- runSim_covs.dat(n, beta, ib, lambda, lambda2, year.entry, MAXT)
    covs.dat    <- tmp$covs.dat
    t.all       <- tmp$t.all
    z2          <- tmp$z2
    year.entry  <- tmp$year.entry
    month.entry <- tmp$month.entry
    c           <- tmp$c
    type        <- tmp$type
    bootflag    <- FALSE
  } else {
    # Take a bootstrap sample
    if (length(boot.rows) != nrow(data_obj$covs.dat)) stop("ERROR 1 with data")
    covs.dat     <- data_obj$covs.dat[boot.rows, , drop=FALSE]
    covs.dat.obs <- data_obj$covs.dat.obs[boot.rows, , drop=FALSE]

    # Let ID values be 1:n to prevent errors later
    covs.dat$ID     <- 1:nrow(covs.dat)
    covs.dat.obs$ID <- 1:nrow(covs.dat.obs)

    year.entry  <- data_obj$year.entry[boot.rows]
    month.entry <- data_obj$month.entry[boot.rows]
    bootflag    <- TRUE
  }

  #simulate mover-stayer model
  for (ims in ims_INDEX){
    if (!bootflag) {
      r_0_sim <-ms.pars[[ims]][1]
      p_00_sim<-ms.pars[[ims]][2]
      p_11_sim<-ms.pars[[ims]][3]
      pi_0_sim<-ms.pars[[ims]][4]
      pi_1_sim<-ms.pars[[ims]][5]
      theta=rnorm(n,mean=0,sd=0.5)
      if (ims>2){
        temp<-log(p_00_sim/(1-p_00_sim))+theta
        p_00_sim<-exp(temp)/(1+exp(temp))
        temp<-log(p_11_sim/(1-p_11_sim))+theta
        p_11_sim<-exp(temp)/(1+exp(temp))
      }
      #Generates simulated observed state of residence using the mover-stayer model
      all_data_sim             <- GenerateAll(r_0_sim, p_00_sim, p_11_sim,n,n.years,pi_0_sim,pi_1_sim)
      sim_data                 <- all_data_sim[[1]]
      stayer_vec_sim_true      <- all_data_sim[[2]]
      rm(all_data_sim)
      gc()

      #True simulated parameters
      init_sim_true <- InitialProbPattern(sim_data,stayer_vec_sim_true)
      tran_sim_true <- TransitionProbPattern(sim_data,stayer_vec_sim_true)
      pi_0_sim_true <- sum(stayer_vec_sim_true == 1)/n
      pi_1_sim_true <- sum(stayer_vec_sim_true == 2)/n

      #removing observations
      sim_data_obs <- IntroduceMissing(sim_data) #Use IntroduceMissing2 if does not allow for 10% missing for the 1st remained data point

      max.month.obs<-apply(sim_data_obs,1,function(x)max(which(!is.na(x))))*12
      max.month.obs[max.month.obs==-Inf]<-NA

      #Simulated observed data
      year.true.calendar<-ceiling((covs.dat$x.true+month.entry)/12)
      temp<-unlist(lapply(1:nrow(sim_data),function(i)sim_data[i,year.true.calendar[i]]==1&covs.dat$delta.true[i]>0))
      VPR.identified<-which(temp==T)
      VPR.unidentified<-which(temp==F)

      if (length(VPR.identified)>0){
        year.identified<-year.true.calendar[VPR.identified]
        sim_data_obs2<-sim_data_obs
        for (i in 1:length(VPR.identified))sim_data_obs2[VPR.identified[i],year.identified[i]]<-1
      } else {
        sim_data_obs2<-sim_data_obs
      }
    
      covs.dat.obs<-covs.dat
      #covs.dat.obs$x.true[VPR.unidentified]<-ceiling(c[VPR.unidentified])
      covs.dat.obs$x.true[VPR.unidentified]<-c[VPR.unidentified]
      covs.dat.obs$delta.true[VPR.unidentified]<-0

      # save sim_data_obs
      # sim_data_obs_list[[ims]] <- sim_data_obs
      sim_data_obs_list[[ims]] <- sim_data_obs2 
    } else {
      sim_data_obs2 <- data_obj$sim_data_obs_list[[ims]]
      # Bootstrap sample
      #sim_data_obs  <- sim_data_obs[boot.rows, , drop=FALSE]
      #sim_data_obs2 <- sim_data_obs 
      sim_data_obs2  <- sim_data_obs2[boot.rows, , drop=FALSE]
    }
    
    start_time<-Sys.time()
    full_data_sim     <- DataSanitation(sim_data_obs2)
    state_pattern_sim <- full_data_sim[[1]]
    freq_vec_sim      <- full_data_sim[[2]]
    rm(full_data_sim)
    gc()

    parameters <- myEM(state_pattern_sim,freq_vec_sim, r_0, p_00, p_11, pi_0, pi_1, epsilon)

    ###########   Generate Masurement Errors ################
    if (!bootflag) {
      type1.1<-type2.2<-rep(NA,n)
      maxt22<-max.month.obs-month.entry
      maxt22[maxt22<0]<-NA
    }
    for (ime1 in ime1_INDEX){
      for (ime2 in 1:1){
        out.id2 <- paste(out.id, ims, ib, ime1, ime2, sep="_") # Output id in files
        if (!bootflag) {
          t2.all=exp(rnorm(n,log(t.all)+z2*me.coef[ime2]/2,me.sd[ime1]))
          type1.1<-ifelse(runif(n)<=sen,1,2)
          type2.2<-ifelse(runif(n)<=spec,2,1)
          type2<-ifelse(type==1,type1.1,type2.2)
          b3<-exp(log(b2)+rnorm(n,me.cen[ime1,ime2],me.sd[ime1]/6))
          c2=(-log(runif(n)))^(1/a2)*b3
          c2<-ifelse(c2<maxt22,c2,maxt22)
          #x2<-ceiling(apply(cbind(t2.all,c2),1,min))
          x2<-apply(cbind(t2.all,c2),1,min)
          delta2<-ifelse(t2.all<c2,type2,0)
    
          covs.dat.obs$x.obs<-x2
          covs.dat.obs$delta.obs<-delta2
          covs.dat$x.obs<-x2
          covs.dat$delta.obs<-delta2
        }
        #####################################################################
        # The following estimators do not depend on the mover-stayer model
        #Gold standard

        event.true<-ifelse(covs.dat$delta.true==1,1,0)
        model1.full=coxph(Surv(x.true,event.true)~.,data=covs.dat[,-c(1:3,5)])

        #VPR data
        event.true<-ifelse(covs.dat.obs$delta.true==1,1,0)
        model1.registry=coxph(Surv(x.true,event.true)~.,data=covs.dat.obs[,-c(1:3,5)])

        #USRT cohort
        event.obs<-ifelse(covs.dat.obs$delta.obs==1,1,0)
        model1.usrt=coxph(Surv(x.obs,event.obs)~.,data=covs.dat.obs[,-c(1,3,4,5)])

        #Hybrid approach
        covs.dat.hybrid<-covs.dat.obs
        covs.dat.hybrid$x.true[covs.dat.obs$delta.true==0&!is.na(covs.dat.obs$delta.obs)]<-covs.dat$x.obs[covs.dat.obs$delta.true==0&!is.na(covs.dat.obs$delta.obs)]
        covs.dat.hybrid$delta.true[covs.dat.obs$delta.true==0&!is.na(covs.dat.obs$delta.obs)]<-covs.dat$delta.obs[covs.dat.obs$delta.true==0&!is.na(covs.dat.obs$delta.obs)]
        event.true<-ifelse(covs.dat.hybrid$delta.true==1,1,0)
        model1.hybrid<-coxph(Surv(x.true,event.true)~.,data=covs.dat.hybrid[,-c(1:3,5)])

        beta.full1<-summary(model1.full)$coef[,1]
        beta.hat1=summary(model1.registry)$coef[,1]
        beta.hat2<-rep(NA,length(beta.hat1))
        #beta.hat2=summary(model1.registry2)$coef[,1]
        #beta.hat.V2=summary(model1.registry.V2)$coef[,1]
        beta.hat.usrt=summary(model1.usrt)$coef[,1]
        beta.hat.hybrid=summary(model1.hybrid)$coef[,1]

        se.beta.full1<-summary(model1.full)$coef[,3]
        se.beta.hat1=summary(model1.registry)$coef[,3]
        se.beta.hat2<-rep(NA,length(se.beta.hat1))
        #se.beta.hat2=summary(model1.registry2)$coef[,3]
        se.beta.hat.usrt=summary(model1.usrt)$coef[,3]
        se.beta.hat.hybrid<-summary(model1.hybrid)$coef[,3]
     
        ####################IDENTIFY VALIDATION SAMPLE and SURVIVAL based on residence status#################################

        #impute missing residence status from the posterior distribution of the mover-stayer model
        p.agree=p.agree2=p.agree3<-rep(NA,M)
        temp.beta.hat.V1<-temp.var.V1<-matrix(NA,nr=M,nc=length(beta.hat1))
        temp.beta.hat.impute1<-temp.beta.hat.impute2<-temp.var.impute1<-temp.var.impute2<-matrix(NA,nr=M*M2,nc=length(beta.hat1))
  
        for (im in 1:M){    
          cat(paste0("im = ", im, "\n"))
          # NOTE: other variables (freq_vec_sim, r_0, p_00, p_11, epsilon) passed in so that the original C code 
          #   written for EM can be re-used. 
          decoded_data <- myDecode(sim_data_obs2, state_pattern_sim, parameters[[1]], parameters[[2]], 
                                   parameters[[3]], parameters[[4]], par.impute,
                                   freq_vec_sim, r_0, p_00, p_11, epsilon)

          end_time<-Sys.time()
          decoded.states<-cbind(1:n,decoded_data)
 
          validation.status   <- unlist(lapply(1:nrow(decoded.states),function(i)decoded.states[i,-1][year.entry[i]+1]==1))
          index.V             <- which(validation.status==T)
          index.nonV          <- which(validation.status==F)
          decoded.states.V    <- decoded.states[index.V,]
          decoded.states.nonV <- decoded.states[index.nonV,]

          if (!bootflag) {
            sim_data.V   <- sim_data[index.V,]
            p.agree[im]  <- mean(decoded.states[,year.entry[1]+2]==sim_data[,year.entry[1]+1])
            p.agree2[im] <- mean(apply(decoded.states.V[,-(1:(year.entry[1]+1))]==sim_data.V[,-(1:year.entry[1])],1,all))
            temp         <- apply(decoded.states.V[,-(1:(year.entry[1]+1))]==sim_data.V[,-(1:year.entry[1])],1,sum)
            p.agree3[im] <- mean(temp>=16) #misclassify no more than 2 residence status
          } else {
            p.agree  <- NA
            p.agree2 <- NA
            p.agree3 <- NA
          }

          covs.dat.V=covs.dat.obs[index.V,]
          covs.dat.nonV=covs.dat.obs[index.nonV,]
  
          ######################################################
          #Censor observation after moving out for cases in the V set
          #covs.decode.dat.V is used for estimate the working model for the predictive distribution

          surv.decode.registry.V=surv.decode3(covs.dat.V$ID,covs.dat.V$x.true,covs.dat.V$delta.true,ceiling(covs.dat.V$x.true/12),
                                              decoded.states.V,1:n.years,year.entry[index.V])
          covs.decode.dat.V=covs.dat.V
          covs.decode.dat.V$x.true=surv.decode.registry.V$surv.time.V
          covs.decode.dat.V$delta.true=surv.decode.registry.V$delta.V

          index.nonV2<-which(validation.status==F&covs.dat.obs$delta.true==0&!is.na(covs.dat.obs$delta.obs))
          covs.dat.nonV2=covs.dat.obs[index.nonV2,]
          futime.nonV2<-covs.dat.nonV2$x.true

          #Generate natural cubic splines of x.obs stratified on delta.obs
          ndf=6
          ns.x.obs=generate.ns(covs.dat$ID,covs.dat.obs$x.obs,covs.dat.obs$delta.obs,ndf)
          #############################################
          #Finding natural splines in the V set for failure time

          subset.V=list()
          uniq.status=sort(unique(covs.decode.dat.V$delta.obs))
          for (i in uniq.status){
            subset.V[[i+1]]=list()

            temp1=covs.decode.dat.V[covs.decode.dat.V$delta.obs==i&!is.na(covs.decode.dat.V$delta.obs),]
            index1=match(ns.x.obs[[i+1]][[1]][,1],temp1[,1])

            for (idf in 1:length(ns.x.obs[[i+1]])){
              ns.mat1=ns.x.obs[[i+1]][[idf]][!is.na(index1),]
              subset.V[[i+1]][[idf]]=data.frame(cbind(temp1[,-c(1:3)],ns.mat1[,-1]))
              temp.nc=ncol(subset.V[[i+1]][[idf]])
              if (idf==1)dimnames(subset.V[[i+1]][[idf]])[[2]][temp.nc]="ns1"
            }
          }

          #Fitting cox model of (x.true,delta.true) on cubin splines of x.obs and covariates stratified by delta.obs

          #spline.fit=coxfit.spline.xobs.cmpr(subset.V)
          #spline.fit=coxfit.spline.xobs.cmpr2(subset.V) #no misclassification of cancer types
          #df.select=spline.fit$df.select
          #df.select<-matrix(3,nr=nrow(df.select),nc=ncol(df.select))
          df.select<-matrix(3,nr=3,nc=2)

          fit.select11<-fit.select12<-fit.select21<-fit.select22<-fit.select31<-fit.select32<-NULL
          dat.V11=subset.V[[1]][[unlist(df.select[1,1])]]
          dat.V11$delta.true<-ifelse(dat.V11$delta.true==1,1,0)
          fit.select11=coxph(Surv(x.true,delta.true)~.,data=dat.V11)

          dat.V12=subset.V[[1]][[unlist(df.select[1,2])]]
          dat.V12$delta.true<-ifelse(dat.V12$delta.true==2,1,0)
          fit.select12=coxph(Surv(x.true,delta.true)~.,data=dat.V12)

          dat.V21=subset.V[[2]][[unlist(df.select[2,1])]]
          dat.V21$delta.true<-ifelse(dat.V21$delta.true==1,1,0)
          fit.select21=coxph(Surv(x.true,delta.true)~.,data=dat.V21)

          if (sen<0.95|spec<0.95){
	    dat.V22=subset.V[[2]][[df.select[2,2]]]
            dat.V22$delta.true<-ifelse(dat.V22$delta.true==2,1,0)
            fit.select22=coxph(Surv(x.true,delta.true)~.,data=dat.V22)
      
            dat.V31=subset.V[[3]][[df.select[3,1]]]
            dat.V31$delta.true<-ifelse(dat.V31$delta.true==1,1,0)
            fit.select31=coxph(Surv(x.true,delta.true)~.,data=dat.V31)
          }

          dat.V32=subset.V[[3]][[unlist(df.select[3,2])]]
          dat.V32$delta.true<-ifelse(dat.V32$delta.true==2,1,0)
          fit.select32=coxph(Surv(x.true,delta.true)~.,data=dat.V32)


          names.fit<-t(matrix(c("fit.select11","fit.select12","fit.select21","fit.select22","fit.select31","fit.select32"),nr=2,nc=3))

          #############################################
          #Validation cohort identified by the moving stayer model

          event.true<-ifelse(covs.decode.dat.V$delta.true==1,1,0) 
          model1.registry.V=coxph(Surv(x.true,event.true)~.,data=covs.decode.dat.V[,-c(1:3,5)])

          temp.beta.hat.V1[im,]=summary(model1.registry.V)$coef[,1]
          temp.var.V1[im,]<-summary(model1.registry.V)$coef[,3]^2
          #beta.hat.V2=summary(model1.registry.V2)$coef[,1]
          #se.beta.hat.V1=summary(model1.registry.V)$coef[,3]
          #se.beta.hat.V2=summary(model1.registry.V2)$coef[,3]

          ################################################################################
          #METHOD 1:

          #Fit kernel-smoothed product-limit surival function for Y.true conditional on Y.obs using V set

          fit.event=fit.cen=list()
          for (i in uniq.status){
            subset=covs.decode.dat.V[covs.decode.dat.V$delta.obs==i&!is.na(covs.decode.dat.V$delta.obs),]
            fit.event[[i+1]]=prodlim(Hist(x.true,delta.true)~x.obs,data=subset) 
          }

          uniq.ftimes<-sort(unique(covs.dat.V$x.true[covs.dat.V$delta.true>0]))
          impute.dat1=out.model1=covs.dat.imputed1=list()

          impute.dat.V<-list() #impute observations in the V set which were censored after moving out of the catchment area
          decoded.states.V<-decoded.states[index.V,]
          year.entry.nonV2<-year.entry[index.nonV2]

	  for (im2 in 1:M2){
            impute.dat1[[im2]]=impute.method.1.cmpr2(par1,covs.dat.nonV2[,c("x.obs","delta.obs","x.true")],fit.event,uniq.ftimes,decoded.states	    [index.nonV2,],year.entry.nonV2)
            covs.dat.imputed1[[im2]]=covs.dat.obs
            covs.dat.imputed1[[im2]][index.nonV2,4]=impute.dat1[[im2]][,1]
            covs.dat.imputed1[[im2]][index.nonV2,5]=impute.dat1[[im2]][,2]  
            index.impute.V<-which(covs.dat.V$delta.true==0&covs.decode.dat.V$x.true!=covs.dat.V$x.true&!is.na(covs.dat.V$delta.obs))
            year.entry.V<-year.entry[index.V]
            impute.dat.V[[im2]]<-impute.method.V.cmpr2(par1,covs.dat.V[index.impute.V,c("x.obs","delta.obs")],covs.decode.dat.V$x.true[index.impute.V],
                                      covs.dat.V$x.true[index.impute.V],fit.event,uniq.ftimes,decoded.states.V[index.impute.V,],year.entry.V[index.impute.V])
            covs.dat.imputed1[[im2]][index.V[index.impute.V],4]=impute.dat.V[[im2]][,1]
            covs.dat.imputed1[[im2]][index.V[index.impute.V],5]=impute.dat.V[[im2]][,2]
            covs.dat.imputed1[[im2]]$delta.true<-ifelse(covs.dat.imputed1[[im2]]$delta.true==1,1,0)
            model1.imputed=coxph(Surv(x.true,delta.true)~.,data=covs.dat.imputed1[[im2]][,-c(1:3)])
            out.model1[[im2]]=summary(model1.imputed)
          }

          if (M2==1){
            temp.beta.hat.impute1[im,]<-out.model1[[1]]$coef[,1] 
            temp.var.impute1[im,]<-out.model1[[1]]$coef[,3]^2
          } else {
            out.coef1=lapply(1:M2, function(i)out.model1[[i]]$coef)
            temp.beta.hat.impute1[((im-1)*M2+1):(im*M2),]=t(do.call(cbind,out.coef1)[,seq(1,M2*5,5)])
            temp.var.impute1[((im-1)*M2+1):(im*M2),]=t(do.call(cbind,out.coef1)[,seq(3,M2*5,5)]^2)
          }
          #var1.1=apply(do.call(cbind,out.coef1)[,seq(3,M*5,5)]^2,1,mean)
          #var2.1=(1+1/M)*apply(do.call(cbind,out.coef1)[,seq(1,M*5,5)],1,var)
          #se.beta.hat.impute1=sqrt(var1.1+var2.1)

          #####################################################################################
          #METHOD 2:

          #Impute Y.true from semi-parametric cox model including B-SPLINE of x.obs and covariates

          #natural cubic spline of observed failure time measured with error stratified by event status
          #degrees of freedom ranges from 1 to 6

          subset.nonV2=list()
          for (i in uniq.status){
            subset.nonV2[[i+1]]=list()
            temp2=covs.dat.nonV2[covs.dat.nonV2$delta.obs==i&!is.na(covs.dat.nonV2$delta.obs),]
            index2=match(ns.x.obs[[i+1]][[1]][,1],temp2[,1])
  
            for (idf in 1:length(ns.x.obs[[i+1]])){
              ns.mat2=ns.x.obs[[i+1]][[idf]][!is.na(index2),]
              subset.nonV2[[i+1]][[idf]]=data.frame(cbind(temp2,ns.mat2[,-1])) #keep x.obs and delta.obs
              temp.nc=ncol(subset.nonV2[[i+1]][[idf]])
              if (idf==1)dimnames(subset.nonV2[[i+1]][[idf]])[[2]][temp.nc]="ns1"
            }
          }

          ns.dat.nonV2=ns.dat.V=list()
          for (i in 1:length(uniq.status)){
            ns.dat.nonV2[[i]]<-list()
            ns.dat.V[[i]]<-list()
            for (j in 1:(length(uniq.status)-1)){
              ns.dat.nonV2[[i]][[j]]<-subset.nonV2[[i]][[df.select[i,j]]]
              ns.dat.V[[i]][[j]]<-subset.V[[i]][[df.select[i,j]]]
            }
          }    

          impute.dat2=out.model2=covs.dat.imputed2=list()
          impute.dat.V2<-list()
          index.impute.V<-which(covs.decode.dat.V$x.true!=covs.dat.V$x.true&covs.dat.V$delta.true==0&!is.na(covs.dat.V$delta.obs)) 

          par22 <- par2
          for (im2 in 1:M2){
            impute.dat2[[im2]]=matrix(NA,nrow(covs.dat.nonV2),2)
            #impute.dat.V2[[im2]]<-matrix(NA,length(index.impute.V),2)
            impute.dat.V2[[im2]]<-matrix(NA,nrow(covs.dat.V),2)
 
            for (i in uniq.status){
              if (is.null(par2)) par22 <- i  # Assuming the index i to be passed in
        
	      index=which(covs.dat.nonV2$delta.obs==i)
              subset.futime<-covs.dat.nonV2$x.true[index]
              #maxt2.index<-maxt2[index.nonV2][index]
              impute.dat2[[im2]][index,]=sample.coxph2.cmpr2(par22,get(names.fit[i+1,1]),get(names.fit[i+1,2]),ns.dat.nonV2[[i+1]][[1]],ns.dat.nonV2[[i+1]][[2]],
              subset.futime,decoded.states[index.nonV2,][index,],year.entry.nonV2[index])
            }
            covs.dat.imputed2[[im2]]=covs.dat.obs
            covs.dat.imputed2[[im2]][index.nonV2,4:5]=impute.dat2[[im2]]
  
            par22 <- par2
            for (i in uniq.status){
              #index.impute.V12<-which(covs.dat.V2$delta.obs==i)
              temp.covs.dat.V<-covs.dat.V[covs.dat.V$delta.obs==i&!is.na(covs.dat.V$delta.obs),]
              temp.covs.decode.dat.V<-covs.decode.dat.V[covs.decode.dat.V$delta.obs==i&!is.na(covs.dat.V$delta.obs),]
              index.impute.V21<-which((covs.dat.V$x.true!=covs.decode.dat.V$x.true&covs.dat.V$delta.true==0)&(covs.dat.V$delta.obs==i&!is.na(covs.dat.V$delta.obs))) #index based on v set
              subset.futime1<-covs.decode.dat.V$x.true[index.impute.V21]
              subset.futime2<-covs.dat.V$x.true[index.impute.V21]
              #maxt2.index<-maxt2[index.V][index.impute.V21]
              index.impute.V22<-which(temp.covs.dat.V$x.true!=temp.covs.decode.dat.V$x.true&temp.covs.dat.V$delta.true==0) #index based on subset of V set with delta.obs=1
              if(length(index.impute.V21)>0){
                if (is.null(par2)) par22 <- i
                impute.dat.V2[[im2]][index.impute.V21,]=sample.coxph.V2.cmpr2(par22,get(names.fit[i+1,1]),get(names.fit[i+1,2]),ns.dat.V[[i+1]][[1]][index.impute.V22,],
                                ns.dat.V[[i+1]][[2]][index.impute.V22,],subset.futime1,subset.futime2,
                                decoded.states.V[index.impute.V21,],year.entry.V[index.impute.V21])
    
              }
            }
            covs.dat.imputed2[[im2]][index.V,4:5][index.impute.V,]<-impute.dat.V2[[im2]][index.impute.V,]
                                                                   
            covs.dat.imputed2[[im2]]$delta.true<-ifelse(covs.dat.imputed2[[im2]]$delta.true==1,1,0)
  
            model2.imputed=coxph(Surv(x.true,delta.true)~.,data=covs.dat.imputed2[[im2]][,-c(1:3)])
            out.model2[[im2]]=summary(model2.imputed)
          }

          if (M2==1) {
            temp.beta.hat.impute2[im,]<-out.model2[[1]]$coef[,1] 
            temp.var.impute2[im,]<-out.model2[[1]]$coef[,3]^2
          } else {
            out.coef2=lapply(1:M2, function(i)out.model2[[i]]$coef)
            temp.beta.hat.impute2[((im-1)*M2+1):(im*M2),]=t(do.call(cbind,out.coef2)[,seq(1,M2*5,5)])
            temp.var.impute2[((im-1)*M2+1):(im*M2),]=t(do.call(cbind,out.coef2)[,seq(3,M2*5,5)]^2)
          }
          #var1.2=apply(do.call(cbind,out.coef)[,seq(3,M*5,5)]^2,1,mean)
          #var2.2=(1+1/M)*apply(do.call(cbind,out.coef)[,seq(1,M*5,5)],1,var)
          #se.beta.hat.impute2=sqrt(var1.2+var2.2)

        } #im=1:M
    
  #      tmp1<-mean(covs.dat$delta.true==0)
  #      tmp2<-mean(covs.dat.obs$delta.true==0)
  #      tmp3<-mean(covs.dat.obs$delta.obs==0,na.rm=T)
  #      tau1<-mean(p.agree)
  #      tau2<-mean(p.agree3)
  #      tmp4<-mean(p.agree2)
  #      pi0<-parameters[[3]]
  #      pi1<-parameters[[4]]
  #      p1<-parameters[[1]][2]
  #      p00<-diag(parameters[[2]])[1]
  #      p11<-diag(parameters[[2]])[2]
 #       out.par<-c(out.id2,ims,me.sd[ime1],me.coef[ime2],tmp1,tmp2,tmp3,tau1,tmp4,tau2,pi0,pi1,p1,p00,p11)
  #     out.par<-c(out.id2,ims,me.sd[ime1],me.coef[ime2],mean(covs.dat$delta.true==0),mean(covs.dat.obs$delta.true==0),mean(covs.dat.obs$delta.obs==0,na.rm=T),
  #                mean(p.agree),mean(p.agree2),mean(p.agree3),parameters[[3]],parameters[[4]],parameters[[1]],diag(parameters[[2]]))
        out.par<-c(out.id2,ims,me.sd[ime1],me.coef[ime2], mean(p.agree),mean(p.agree3),parameters[[3]],parameters[[4]],parameters[[1]],diag(parameters[[2]]))
        write(out.par,file=out.files[3],sep="\t",ncol=length(out.par),append=T)

        ################################################################################
        #Calculate multiple-imputation based estimate of beta's and SEs by Rubin's rule
        beta.hat.V1<-apply(temp.beta.hat.V1,2,mean)
        if (M==1) se.beta.hat.V1<-sqrt(temp.var.V1) else { 
          var1.2<-apply(temp.var.V1,2,mean)
          var2.2<-(1+1/M)*apply(temp.beta.hat.V1,2,var)
          se.beta.hat.V1<-sqrt(var1.2+var2.2)
        }

        beta.hat.impute1<-apply(temp.beta.hat.impute1,2,mean)
        var1.2<-apply(temp.var.impute1,2,mean)
        var2.2<-(1+1/(M*M2))*apply(temp.beta.hat.impute1,2,var)
        se.beta.hat.impute1<-sqrt(var1.2+var2.2)

        beta.hat.impute2<-apply(temp.beta.hat.impute2,2,mean)
        var1.2<-apply(temp.var.impute2,2,mean)
        var2.2<-(1+1/(M*M2))*apply(temp.beta.hat.impute2,2,var)
        se.beta.hat.impute2<-sqrt(var1.2+var2.2)

        ####Calculate coverage probability ####
        BETA         <- beta[ib,]
        covp.full1   <- getCoverageProb(BETA, beta.full1,       se.beta.full1,       z=1.96)
        covp.hat1    <- getCoverageProb(BETA, beta.hat1,        se.beta.hat1,        z=1.96)
        covp.hybrid  <- getCoverageProb(BETA, beta.hat.hybrid,  se.beta.hat.hybrid,  z=1.96)
        covp.usrt    <- getCoverageProb(BETA, beta.hat.usrt,    se.beta.hat.usrt,    z=1.96)
        covp.V1      <- getCoverageProb(BETA, beta.hat.V1,      se.beta.hat.V1,      z=1.96)
        covp.impute1 <- getCoverageProb(BETA, beta.hat.impute1, se.beta.hat.impute1, z=1.96)
        covp.impute2 <- getCoverageProb(BETA, beta.hat.impute2, se.beta.hat.impute2, z=1.96)

        out.est=c(out.id2,ims,NA,beta[ib,],me.sd[ime1],me.coef[ime2],beta.full1,beta.hat1,beta.hat.hybrid,beta.hat.usrt,
                  beta.hat.V1,NA,NA,beta.hat.impute1,NA,NA,beta.hat.impute2,NA,NA)
        #          beta.hat.em)
        out.se=c(out.id2,ims,NA,beta[ib,],me.sd[ime1],me.coef[ime2],se.beta.full1,se.beta.hat1,se.beta.hat.hybrid,
                 se.beta.hat.usrt,se.beta.hat.V1,NA,NA,se.beta.hat.impute1,NA,NA,se.beta.hat.impute2,NA,NA,covp.full1,covp.hat1,covp.hybrid,covp.usrt,
                 covp.V1,covp.impute1,covp.impute2)
        write(out.est,file=out.files[1],sep="\t",ncol=length(out.est),append=T)
        write(out.se, file=out.files[2],sep="\t",ncol=length(out.se),append=T)

      }#ime1
    }#ime2
  }#ims 1:4

  ret <- list(covs.dat=covs.dat, covs.dat.obs=covs.dat.obs, year.entry=year.entry, month.entry=month.entry,
              sim_data_obs_list=sim_data_obs_list)
  ret

} # END: runSim_main

getCoverageProb <- function(beta, betahat, sehat, z=1.96) {
 
  # beta, betahat, and sehat are vectors

  low  <- betahat - z*sehat
  upp  <- betahat + z*sehat
  covp <- ifelse(low < beta & upp > beta, 1, 0)

  covp

}

getCoverageProb_perc <- function(beta, betahat, perc=c(0.025, 0.975)) {
 
  # beta is a vector
  # betahat is a matrix

  n    <- ncol(betahat)
  covp <- rep(NA, n)
  for (i in 1:n) {
    qq      <- quantile(betahat[, i], probs=perc, na.rm=TRUE)
    covp[i] <- ifelse((qq[1] < beta[i]) && (qq[2] > beta[i]), 1, 0)
  }
  covp

}

getBootSE_file <- function(obsFile, bootFiles, parms, outFile, perc=c(0.025, 0.975)) {

  # Column names for the betahats
  nc.beta  <- ncol(parms$beta)
  tmp      <- c("beta.GS", "beta.IR", "beta.Hybrid", "beta.SR", 
                "beta.Step_1", "beta.NP", "beta.Proposed")
  betavars <- paste(rep(tmp, each=nc.beta), rep(1:nc.beta, times=length(tmp)), sep=".") 
  ib       <- parms$ib
  BETA     <- parms$beta[ib, , drop=TRUE]
  nbootf   <- length(bootFiles)
  
  # Read in the file containing the betahats
  x0  <- read.table(obsFile, header=1, sep="\t", stringsAsFactors=FALSE)
  x0  <- unique(x0)
  if (!any(grepl("boot0", x0[, 1], fixed=TRUE))) stop("ERROR with obsFile")
  if (!all(betavars %in% colnames(x0))) stop("ERROR with obsFile")
  x0[, 1]      <- gsub("_boot0", "", x0[, 1, drop=TRUE], fixed=TRUE)
  ids          <- x0[, 1, drop=TRUE]
  if (any(duplicated(ids))) stop("ERROR with obsFile")
  rownames(x0) <- ids

  # Read in the file containing the bootstrap results
  x <- NULL
  for (i in 1:nbootf) {
    y  <- read.table(bootFiles[i], header=1, sep="\t", stringsAsFactors=FALSE)
    if (!nrow(y)) next
    if (any(grepl("boot0", y[, 1], fixed=TRUE))) stop("ERROR with bootFile")
    tmp <- colnames(y)
    if (!all(betavars %in% tmp)) stop("ERROR with bootFile")
    if (i < 2) {
      x <- y
    } else {
      x <- rbind(x, y)
    }
  }
  if (!length(x)) stop("ERROR: no bootstraps")
  y <- NULL
  gc()
  cx <- colnames(x)

  for (v in betavars) {
    x0[, v] <- as.numeric(x0[, v, drop=TRUE])
    x[, v]  <- as.numeric(x[, v, drop=TRUE])  
  }
  nrx   <- nrow(x)
  nbeta <- length(betavars)/nc.beta
  nids  <- length(ids)

  # Get bootstrap standard errors for each replication, remove "boot1", "boot2", ... etc
  tmp    <- matrix(unlist(strsplit(x[, 1, drop=TRUE], "_", fixed=TRUE)), nrow=nrx, ncol=6, byrow=TRUE)
  idboot <- myPasteCols(tmp[, -2, drop=FALSE], sep="_") 

  # Initialize matrix for results
  col1          <- gsub("_boot", "", cx[1], fixed=TRUE)
  secols        <- paste0("boot.se.", betavars)
  covpcols      <- paste0("covp.", betavars) 
  covpcols2     <- paste0("covpPerc.", betavars) # Percentile method 
  nonmisscols   <- paste0("nboot.", betavars) # number of non-missing
  nacols        <- paste0("NA.", 1:99)
  cols2         <- cx[!(cx %in% c(betavars, nacols, cx[1]))]
  tmp           <- c(col1, cols2, secols, covpcols, covpcols2, nonmisscols)
  ret           <- matrix(data=NA, nrow=nids, ncol=length(tmp)) 
  colnames(ret) <- tmp
  ret[, 1]      <- ids
  ret[, cols2]  <- as.matrix(x0[, cols2]) 

  for (i in 1:nids) {
    id  <- ids[i]
    tmp <- idboot %in% id
    if (!any(tmp)) next

    x2  <- x[tmp, , drop=FALSE]
    b   <- 0
    for (j in 1:nbeta) {
      a                    <- b + 1
      b                    <- a + nc.beta - 1
      ab                   <- a:b
      bvars                <- betavars[ab]
      betahat              <- unlist(x0[id, bvars]) 
      bootse               <- rep(NA, nc.beta)
      nonmiss              <- bootse
      for (k in 1:nc.beta) {
        vec        <- x2[, bvars[k], drop=TRUE]
        bootse[k]  <- sd(vec, na.rm=TRUE) 
        nonmiss[k] <- sum(is.finite(vec))
      }
      ret[i, secols[ab]]      <- bootse 
      ret[i, covpcols[ab]]    <- getCoverageProb(BETA, betahat, bootse, z=1.96)
      ret[i, nonmisscols[ab]] <- nonmiss

      # Perentile method of coverage probs
      ret[i, covpcols2[ab]] <- getCoverageProb_perc(BETA, x2[, bvars, drop=FALSE], perc=perc)
    }
  }

  if (length(outFile)) write.table(ret, file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  ret
}

getBootSE_file2 <- function(obsFile, bootFiles, parms, outFile, perc=c(0.025, 0.975)) {

  # Column names for the parameters in the mover-stayer model
  tmp <-c("tau1","tau2","pi0","pi1","p1","p00","p11")
  betavars <- tmp
  # nc.beta  <- ncol(parms$beta)
  nc.beta<-length(betavars)
 # ib       <- parms$ib
 # BETA     <- parms$beta[ib, , drop=TRUE]
  nbootf   <- length(bootFiles)
  
  # Read in the file containing the betahats
  x0  <- read.table(obsFile, header=1, sep="\t", stringsAsFactors=FALSE)
  x0  <- unique(x0)
  if (!any(grepl("boot0", x0[, 1], fixed=TRUE))) stop("ERROR with obsFile")
  if (!all(betavars %in% colnames(x0))) stop("ERROR with obsFile")
  x0[, 1]      <- gsub("_boot0", "", x0[, 1, drop=TRUE], fixed=TRUE)
  ids          <- x0[, 1, drop=TRUE]
  if (any(duplicated(ids))) stop("ERROR with obsFile")
  rownames(x0) <- ids

  # Read in the file containing the bootstrap results
  x <- NULL
  for (i in 1:nbootf) {
    y  <- read.table(bootFiles[i], header=1, sep="\t", stringsAsFactors=FALSE)
    if (!nrow(y)) next
    if (any(grepl("boot0", y[, 1], fixed=TRUE))) stop("ERROR with bootFile")
    tmp <- colnames(y)
    if (!all(betavars %in% tmp)) stop("ERROR with bootFile")
    if (i < 2) {
      x <- y
    } else {
      x <- rbind(x, y)
    }
  }
  if (!length(x)) stop("ERROR: no bootstraps")
  y <- NULL
  gc()
  cx <- colnames(x)

  for (v in betavars) {
    x0[, v] <- as.numeric(x0[, v, drop=TRUE])
    x[, v]  <- as.numeric(x[, v, drop=TRUE])  
  }
  nrx   <- nrow(x)
  nbeta <- length(betavars)/nc.beta
  nids  <- length(ids)

  # Get bootstrap standard errors for each replication, remove "boot1", "boot2", ... etc
  tmp    <- matrix(unlist(strsplit(x[, 1, drop=TRUE], "_", fixed=TRUE)), nrow=nrx, ncol=6, byrow=TRUE)
  idboot <- myPasteCols(tmp[, -2, drop=FALSE], sep="_") 

  # Initialize matrix for results
  col1          <- gsub("_boot", "", cx[1], fixed=TRUE)
  secols        <- paste0("boot.se.", betavars)
  #covpcols      <- paste0("covp.", betavars) 
  #covpcols2     <- paste0("covpPerc.", betavars) # Percentile method 
  nonmisscols   <- paste0("nboot.", betavars) # number of non-missing
  nacols        <- paste0("NA.", 1:99)
  cols2         <- cx[!(cx %in% c(betavars, nacols, cx[1]))]
  #tmp           <- c(col1, cols2, secols, covpcols, covpcols2, nonmisscols)
  tmp           <- c(col1, cols2, secols, nonmisscols)
  ret           <- matrix(data=NA, nrow=nids, ncol=length(tmp)) 
  colnames(ret) <- tmp
  ret[, 1]      <- ids
  ret[, cols2]  <- as.matrix(x0[, cols2]) 

  for (i in 1:nids) {
    id  <- ids[i]
    tmp <- idboot %in% id
    if (!any(tmp)) next

    x2  <- x[tmp, , drop=FALSE]
    b   <- 0
    for (j in 1:nbeta) {
      a                    <- b + 1
      b                    <- a + nc.beta - 1
      ab                   <- a:b
      bvars                <- betavars[ab]
      betahat              <- unlist(x0[id, bvars]) 
      bootse               <- rep(NA, nc.beta)
      nonmiss              <- bootse
      for (k in 1:nc.beta) {
        vec        <- x2[, bvars[k], drop=TRUE]
        bootse[k]  <- sd(vec, na.rm=TRUE) 
        nonmiss[k] <- sum(is.finite(vec))
      }
      ret[i, secols[ab]]      <- bootse 
      #ret[i, covpcols[ab]]    <- getCoverageProb(BETA, betahat, bootse, z=1.96)
      ret[i, nonmisscols[ab]] <- nonmiss

      # Perentile method of coverage probs
      #ret[i, covpcols2[ab]] <- getCoverageProb_perc(BETA, x2[, bvars, drop=FALSE], perc=perc)
    }
  }

  if (length(outFile)) write.table(ret, file=outFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  ret
}

check_parms <- function(op) {

  # par1 and par2 are used in:
  # impute.method.1.cmpr2(par1,…) 
  # impute.method.V.cmpr2(par1,…)
  # sample.coxph2.cmpr2(par2,…)
  # sample.coxph.V2.cmpr2(par2,…) 

  if (!length(op)) op <- list()
  if (!is.list(op)) stop("ERROR: op must be a list")

  # Valid parameter names
  valid   <- c("M", "M2", "n", "beta", 
               "lambda", "lambda2", 
               "a", "b", "a2", "b2", 
               "MAXT", "me.sd", "me.coef", 
               "me.cen", "ms.pars", "n.years",
               "r_0", "p_00", "p_11", "pi_0", "pi_1",
               "epsilon", "sen", "spec", "ib",
               "par1", "par2", "par.impute",
               "ims_INDEX", "ime1_INDEX")

  # Default parameter values
  default <- list(20, 1, 5000, t(matrix(c(0,0,0.2,-0.2,0.4,-0.4),nr=2,nc=3)), 
                  exp(-9.1), NULL,
                  1, exp(9), 1/0.22, exp(5), 
                  300, c(1,1.5,2), c(0,0.1), 
                  t(matrix(c(-0.05,-0.075,-0.5,-0.5,-1,-1),nr=2,nc=3)), NULL, 25,
                  .3932528, .9712229, .979527, .2201757, .5087876,
                  1e-7, 0.99, 0.99, 1,
                  NULL, NULL, TRUE,
                  NULL, NULL)

  nms <- names(op)
  if (length(nms)) {
    tmp <- !(nms %in% valid)
    if (any(tmp)) {
      err <- paste(nms[tmp], collapse=", ", sep="")
      msg <- paste("ERROR: the options ", err, " are not valid option names")
      stop(msg)
    }
  }

  op <- default.list(op, valid, default)
  if (is.null(op[["lambda2", exact=TRUE]])) op$lambda2 <- 3.5*op$lambda
  
  #Initial parameters for simulation of the mover-stayer model
  #11/9/2020: three sets of initial parameters for simulation of the mover-stayer model
  if (is.null(op[["ms.pars", exact=TRUE]])) {
    ms.pars<-list()
    ms.pars[[1]]<-c(0.39,0.97,0.98,0.2,0.4)
    ms.pars[[2]]<-c(0.39,0.9,0.9,0.2,0.4)
    ms.pars[[3]]<-ms.pars[[1]] #add random effects to the transition probability
    ms.pars[[4]]<-ms.pars[[2]] #add random effects to the transition probability
    op$ms.pars <- ms.pars
  }

  op

} # END: check_parms






#define IARG_DEBUG 0
#define IARG_data_nr 1
#define IARG_data_nc 2
#define IARG_freqVecLen 3
#define IARG_maxiter 4
#define IARG_print 5
#define IARG_kLen 6

#define DARG_r_0 0
#define DARG_p_00 1
#define DARG_p_11 2
#define DARG_pi_0 3
#define DARG_pi_1 4
#define DARG_epsilon 5
#define DARG_machine_eps 6

get_args <- function(data, freq_vec, r_0, p_00, p_11, pi_0, pi_1, epsilon) {
 
  reti <- c(0, nrow(data), ncol(data), length(freq_vec), 100000, 1, 2)
  if (length(reti) != 7) stop("ERROR1")

  retd <- c(r_0, p_00, p_11, pi_0, pi_1, epsilon, .Machine$double.eps)
  if (length(retd) != 7) stop("ERROR2")

  list(iargs=reti, dargs=retd)

}

myEM <- function(data, freq_vec, r_0, p_00, p_11, pi_0, pi_1, epsilon) {

  tmp       <- get_args(data, freq_vec, r_0, p_00, p_11, pi_0, pi_1, epsilon)
  iargs     <- tmp$iargs
  dargs     <- tmp$dargs
  kvec      <- numOfStates(data)
  klen      <- length(kvec)
  ret_code  <- as.integer(-1)
  ret_init  <- as.numeric(rep(-1, klen))
  ret_trans <- as.numeric(rep(-1, klen*klen))
  ret_pi0   <- as.numeric(-1)
  ret_pi1   <- as.numeric(-1)

  MISSING   <- 9
  datavec   <- as.vector(t(data))
  tmp       <- !is.finite(datavec)
  if (any(tmp)) datavec[tmp] <- MISSING

  tmp <- .C("MY_C_EM", as.integer(iargs), as.numeric(dargs), 
                    as.integer(datavec), as.integer(freq_vec), as.integer(kvec),
                    ret_code=ret_code, ret_init=ret_init, ret_trans=ret_trans,
                    ret_pi0=ret_pi0, ret_pi1=ret_pi1)
  
  ret_code  = tmp$ret_code
  if (ret_code != 1) stop("EM algorithm did not converge")
  ret_init  = tmp$ret_init
  ret_trans = matrix(data=tmp$ret_trans, nrow=klen, ncol=klen, byrow=TRUE)
  ret_pi0   = tmp$ret_pi0
  ret_pi1   = tmp$ret_pi1

  list(ret_init, ret_trans, ret_pi0, ret_pi1)
}


my_newt <- function(temp.t.time, temp.t.surv, rv1) {

  # newt=unlist(lapply(1:ncol(temp.t.surv),function(i)temp.t.time[which(temp.t.surv[,i]<=rv1[i])[1]]))

  nc <- ncol(temp.t.surv)
  if (nc < 2) {
    ret <- temp.t.time[which(temp.t.surv[, 1] <= rv1[1])[1]]
  } else {
    nr  <- nrow(temp.t.surv)
    z   <- temp.t.surv <= matrix(rep(rv1, times=nr), nrow=nr, ncol=nc, byrow=TRUE)
    z[is.na(z)] <- FALSE
    v   <- max.col(t(z), ties.method="first")
    ii  <- cbind(v, 1:nc)
    ret <- temp.t.time[v]
    tmp <- !z[ii]
    if (any(tmp)) ret[tmp] <- NA 
  }

  ret
}

my_captured.index <- function(impute.year, year0.index, decode.index, impute.status) {

  #captured.index<-unlist(lapply(1:nrow(decode.index),function(i)ifelse(decode.index[i,-1][impute.year[i]+year0.index[i]]==1&impute.status[i]!=0,1,0)))

  nr <- nrow(decode.index)
  if (nr < 2) {
    ret <- ifelse(decode.index[1,-1][impute.year[1]+year0.index[1]]==1&impute.status[1]!=0,1,0)
  } else {
    tmp <- impute.year + year0.index
    tmp <- cbind(1:nr, tmp)
    mat <- decode.index[, -1, drop=FALSE]
    tmp <- (mat[tmp] == 1) & (impute.status != 0)
    tmp[is.na(tmp)] <- FALSE
    ret <- as.numeric(tmp)
  }
  ret
}

my_surv.t <- function(temp.t.time, temp.t.surv, futime12) {

  # surv.t<-unlist(lapply(1:length(futime12),function(i)temp.t.surv[max(which(temp.t.time<=futime12[i])),i]))

  len <- length(futime12) 
  if (len < 2) {
    ret <- temp.t.surv[max(which(temp.t.time<=futime12[1])),1]
  } else {
    ntime       <- length(temp.t.time)
    z           <- temp.t.time <= matrix(data=rep(futime12, each=ntime), nrow=ntime, ncol=len, byrow=FALSE)
    z[is.na(z)] <- FALSE
    v           <- max.col(t(z), ties.method="last")
    ii          <- cbind(v, 1:len)
    ret         <- temp.t.surv[ii]
    tmp         <- !z[ii]
    if (any(tmp)) ret[tmp] <- NA 
  }

  ret
}

my_temp.inc <- function(temp.cinc) {

  #for (j in 1:length(temp.cinc)){
  #if (is.list(temp.cinc[[j]])==T)
  #temp.inc[[j]]<-do.call(rbind,lapply(1:length(temp.cinc[[j]]),function(k)temp.cinc[[j]][[k]]-c(0,temp.cinc[[j]][[k]][1:(length(temp.cinc[[j]][[k]])-1)]))) 
  #else temp.inc[[j]]<-matrix(temp.cinc[[j]]-c(0,temp.cinc[[j]][1:(length(temp.cinc[[j]])-1)]),nr=1)
  #}

  temp.inc <- list()
  
  for (j in 1:length(temp.cinc)) {
    tmp <- temp.cinc[[j]]
    nr  <- length(tmp)
    if (is.list(tmp)) {
      len <- length(tmp[[1]])
      mat <- matrix(data=flatten_dbl(tmp), nrow=nr, ncol=len, byrow=TRUE)
    } else {
      len <- length(tmp)
      mat <- matrix(data=tmp, nrow=1, ncol=len, byrow=TRUE)
    }  
    cols2         <- 2:len 
    tmp           <- mat
    tmp[, cols2]  <- mat[, cols2, drop=FALSE] - mat[, cols2-1, drop=FALSE]
    temp.inc[[j]] <- tmp  
  }
  temp.inc
}

my_temp.inc0 <- function(temp.cinc) {

  temp.inc <- list()

  for (j in 1:length(temp.cinc)) {
    tmp <- temp.cinc[[j]] 
    len <- length(tmp)
    if (is.list(tmp)) {
      nc   <- length(tmp[[1]])
      mat  <- matrix(data=NA, nrow=len, ncol=nc)
      cols <- 1:(nc - 1)
      for (k in 1:len) {
        vec      <- tmp[[k]] 
        mat[k, ] <- vec - c(0, vec[cols]) 
      }
      temp.inc[[j]] <- mat
    } else {
      temp.inc[[j]]<-matrix(tmp - c(0, tmp[1:(len-1)]), nrow=1)
    }
  }
  
  temp.inc
}

default.list <- function(inList, names, default, error=NULL,
                         checkList=NULL) {

  # inList      List
  # names       Vector of names of items in inList
  # default     List of default values to assign if a name is not found
  #             The order of default must be the same as in names.
  # error       Vector of TRUE/FALSE if it is an error not to have the
  #             name in the list. 
  #             The default is NULL
  # checkList   List of valid values for each name.
  #             Use NA to skip a list element.
  #             The default is NULL

  n1 <- length(names)
  n2 <- length(default)
  if (n1 != n2) stop("ERROR: in calling default.list")

  if (is.null(error)) {
    error <- rep(0, times=n1)
  } else if (n1 != length(error)) {
    stop("ERROR: in calling default.list")
  }

  if (!is.null(checkList)) {
    if (n1 != length(checkList)) stop("ERROR: in calling default.list")
    checkFlag <- 1
  } else {
    checkFlag <- 0
  } 

  if (is.null(inList)) inList <- list()

  listNames <- names(inList)
  for (i in 1:n1) {
    if (!(names[i] %in% listNames)) {
      if (!error[i]) {
        inList[[names[i]]] <- default[[i]]
      } else {
        temp <- paste("ERROR: the name ", names[i], " was not found", sep="")
        stop(temp)
      }
    } else if (checkFlag) {
      temp <- checkList[[i]]
      if (!all(is.na(temp))) {
        if (!all(inList[[names[i]]] %in% checkList[[i]])) {
          temp <- paste("ERROR: the name '", names[i], 
                      "' has an invalid value", sep="")
          stop(temp)
        }
      }
    }
  }

  inList

} # END: default.list

getSwarmFiles <- function(cfiles, call.dir, filePrefix, nswarm=1, BIOWULF=1) {

  # appcall files
  appfiles <- gsub(".R", "", cfiles, fixed=TRUE)
  lfiles   <- paste(cfiles, ".log", sep="")
  for (i in 1:length(appfiles)) {
      f   <- appfiles[i]
      fid <- file(f, "w")
      cat("#!/bin/bash \n", file=fid)
      #cat("module load gcc/4.8.4\n", file=fid)
      cat("module load R/3.6 \n", file=fid)
      #cat("module load plink \n", file=fid)
      str <-  paste("R --vanilla <  ", f, ".R > ", lfiles[i], sep="")
      write(str, file=fid)
      close(fid)
  }

  out <- paste(call.dir, "Rjobs_", filePrefix, sep="")
  if (BIOWULF) {
    str <- paste("bash ", appfiles, sep="")
  } else {
    str <- paste("qsub ", qsub.op, " ", appfiles, sep="")
  }
  N   <- length(str)
  b   <- 0
  M   <- ceiling(N/nswarm)
  for (i in 1:nswarm) {
    out  <- paste(call.dir, "Rjobs", i, "_", filePrefix, sep="")
    a    <- b + 1
    b    <- a + M - 1
    if (b > N) b <- N
    if (a > N) break
    str2 <- c("#!/bin/bash", str[a:b])
    write(str2, file=out)
  }

  NULL

}

getStartBootIndex <- function(f) {

  if (!file.exists(f)) return(1)
  x  <- scan(f, what="char", sep="\n", quiet=TRUE)
  nx <- length(x)
  if (nx < 2) return(1)
  x   <- x[nx]
  x   <- unlist(strsplit(x, "\t", fixed=TRUE))
  vec <- unlist(strsplit(x, "_", fixed=TRUE))
  str <- vec[2]
  str <- gsub("boot", "", str, fixed=TRUE)
  ret <- as.numeric(str) + 1
  if (!is.finite(ret)) stop("ERROR")
  ret

}

runSims_genSwarmFiles <- function(jobfiles.dir0, out.dir0, sourceFiles, parms, seeds.rep, seeds.boot=1:500,
                                  nRepPerJob=1, nBootPerJob=NULL, op=NULL, filePrefix="", nswarm=0, append=0) {
  
  # jobfiles.dir        Folder to create the swarm job files to submit
  # out.str             Output string that includes the complete path. Ex: /data/wheelerwi/simulation/out/sim1
  # sourceFiles         Character vector of source files
  # parms               List of parameters
  # seeds.rep           Numeric vector of replication seeds. Each replication should have its own seed.  
  # seeds.boot          Numeric vector of bootstrap seeds. The default is 1:500, ie, 500 bootstraps will be done.
  # nRepPerJob          Number of replications per job file. The default is 1.
  # nBootPerJob         Number of bootstraps per job file. Use 0 for no bootstrapping. The default is length(seeds.boot).
  # filePrefix          File prefix for the job files. The default is "".

  ims_INDICES  <- 1:length(parms$ms.pars)
  #ime1_INDICES <- 1:3
  ime1_INDICES<-1:length(parms$me.sd)
  models       <- getSeqsFromList(list(ms=ims_INDICES, ime1=ime1_INDICES))
print(models)
  nmodels      <- nrow(models)
  if (is.null(op)) op <- list()

  nrep  <- length(seeds.rep)
  if (!nrep) stop("ERROR with seeds.rep")
  if (nRepPerJob < 1) stop("ERROR with nRepPerJob")
  nboot <- length(seeds.boot)
  if (nboot) {
    if (is.null(nBootPerJob)) nBootPerJob <- nboot
    if (nBootPerJob > nboot) nBootPerJob <- nboot
    if (nBootPerJob < 1) stop("ERROR with nBootPerJob")
    if (nboot == nBootPerJob) {
      op$boot.summary <- 1
    } else {
      op$boot.summary <- 0
    }
  }
  
  len <- nchar(jobfiles.dir0)
  if (substr(jobfiles.dir0, len, len) != "/") jobfiles.dir0 <- paste0(jobfiles.dir0, "/")
  len <- nchar(out.dir0)
  if (substr(out.dir0, len, len) != "/") out.dir0 <- paste0(out.dir0, "/")
  if (!file.exists(jobfiles.dir0)) system(paste0("mkdir ", jobfiles.dir0))
  if (!file.exists(out.dir0)) system(paste0("mkdir ", out.dir0))
  cdir0  <- jobfiles.dir0
  CFILES <- NULL
  for (model in 1:nmodels) {
    ims_INDEX  <- models[model, 1]
    ime1_INDEX <- models[model, 2] 
    id0        <- paste0("ms.", ims_INDEX, "_ime1.", ime1_INDEX)
    cdir       <- paste0(jobfiles.dir0, "call_", id0, "/")
    odir       <- paste0(out.dir0, "out_", id0, "/")
    if (!file.exists(cdir)) system(paste0("mkdir ", cdir))
    if (!file.exists(odir)) system(paste0("mkdir ", odir))
    parms$ims_INDEX  <- ims_INDEX
    parms$ime1_INDEX <- ime1_INDEX

    cfiles     <- NULL
    a.rep      <- 0
    b.rep      <- 0
    a.boot     <- 0
    b.boot     <- 0
    i.rep      <- 0
    i.boot     <- 0
    jobno      <- 0
    while (1) {
      jobno  <- jobno + 1 
      a.rep  <- b.rep + 1
      b.rep  <- a.rep + nRepPerJob - 1
      if (b.rep > nrep) b.rep <- nrep
      if (a.rep > nrep) {
        break
      }

      if (nboot) {
        a.boot  <- b.boot + 1
        b.boot  <- a.boot + nBootPerJob - 1
        if (b.boot > nboot) b.boot <- nboot
        if (a.boot > nboot) {
          # Re-initialize
          a.boot <- 1
          b.boot <- a.boot + nBootPerJob - 1
        }
      }

      id        <- paste0("job", jobno)
      out.str2  <- paste0(odir, id)
      outF      <- paste0(out.str2, "_DUMY5_boot")
      op$append <- 0
      op$out.id.boot <- 1
      if (append) {
        if (file.exists(outF)) next
        # Make sure files exist 
        ff <- paste0(out.str2, "_DUMY", 1:3)
        if (all(file.exists(ff))) {
          # Get starting index
          a.boot <- getStartBootIndex(paste0(out.str2, "_DUMY1_boot"))
          if (a.boot > nboot) stop("ERROR1 a.boot")
          if (a.boot > b.boot) stop("ERROR2 a.boot")
          if (a.boot > 1) {
            op$append      <- 1
            op$out.id.boot <- a.boot
          }
        }
      } 
      temp     <- paste(cdir, filePrefix, id, ".R", sep="")
      fid      <- file(temp, "w")
      cfiles   <- c(cfiles, temp)

      # replication seeds
      rseeds <- seeds.rep[a.rep:b.rep] 

      # bootstrap seeds
      if (nboot) {
        bseeds <- seeds.boot[a.boot:b.boot] 
      } else {
        bseeds <- NULL
      }

      # Source files
      for (f in sourceFiles) {
        f    <- getLoadSourceFileStr(f)
        temp <- paste(f, '\n', sep="")
        cat(temp, file=fid)
      }
   
      writeToFile(out.str2, "out.str2", fid)
      writeToFile(rseeds, "rseeds", fid)
      writeToFile(bseeds, "bseeds", fid)
      writeToFile(parms, "parms", fid)
      writeToFile(op, "op", fid)

      temp <- "runSims(out.str2, rseeds, parms, seeds.boot=bseeds, op=op)\n"
      cat(temp, file=fid)
      close(fid)
    }
    if (!nswarm) {
      getSwarmFiles(cfiles, cdir, filePrefix, nswarm=1)
    } else {
      CFILES <- c(CFILES, cfiles)
    }
  }
  if (nswarm) getSwarmFiles(CFILES, cdir0, filePrefix, nswarm=nswarm)
   
  NULL

} 

getLoadSourceFileStr <- function(f) {

  len <- nchar(f)
  lf  <- tolower(f)
  str <- substr(lf, len-1, len)
  if (str == ".r") {
    ret <- paste0('source("', f, '")')
  } else if (str == "so") {
    ret <- paste0('dyn.load("', f, '")')
  } else {
    stop(paste0("ERROR with source file ", f))
  }
  ret

}

writeToFile <- function(obj, name, fid) {

  if (is.null(obj)) {
    temp <- paste(name, " <- NULL \n", sep="")
    cat(temp, file=fid)
  } else if (isString(obj)) {
    temp <- paste(name, ' <- "', obj, '"\n', sep="")
    cat(temp, file=fid)
  } else if (is.matrix(obj)) {
    nr   <- nrow(obj)
    nc   <- ncol(obj)
    vec  <- as.vector(obj)
    str  <- paste0("c(", paste0(vec, collapse=","), ")")
    temp <- paste0(name, " <- matrix(", str, ",nrow=", nr, ",ncol=", nc, ")\n")
    cat(temp, file=fid)
  } else {
    if (!is.list(obj)) {
      if (length(obj) == 1) {
        temp <- paste(name, " <- ", obj, "\n", sep="")
        cat(temp, file=fid)
      } else {
        genfile.vec(obj, name, fid)
      }
    } else {
      genfile.list(obj, name, fid)
    }
  }
  NULL
}

isString <- function(obj) {
  ret <- (length(obj) == 1) && is.character(obj)
  ret
}

# Function to define list in swarm generator files
genfile.list <- function(inList, listName, fid) {

  # If the field is inList is a function or family,
  # put it in quotes and set the comment to "FUNCTION".

  if (is.null(inList)) {
    temp <- paste(listName, " <- NULL\n", sep="")
    cat(temp, file=fid)
    return(NULL)
  } else if (!length(inList)) {
    temp <- paste(listName, " <- list()\n", sep="")
    cat(temp, file=fid)
    return(NULL)
  }

  names <- names(inList)
  llen  <- length(inList)
  if (is.null(names)) {
    flag  <- 1   
    names <- 1:llen
    str1  <- '[['
    str2  <- ']] <- '
  } else {
    flag <- 0
    str1 <- '[["'
    str2 <- '"]] <- '
  }
  temp <- nchar(names) < 1
  m    <- sum(temp) 
  if (m == llen) {
    names <- paste("x", 1:llen, sep="")
  } else if (m) {
    new <- getUniqueVarName(names[!temp], alen=2, nlen=2)
    names[temp] <- paste(new, "_", 1:m, sep="")
  }

  temp  <- paste("\n ", listName, " <- list() \n", sep="")
  cat(temp, file=fid)

  if (!llen) return(NULL)
  for (name in names) {
    temp <- inList[[name]]
    cmm  <- comment(temp)
    if (is.null(cmm)) cmm <- "NULL"

    # For a list
    if (is.list(temp)) {
      if (flag) {
        name2 <- paste("list", name, sep="")
      } else {
        name2 <- name
      }
      genfile.list(temp, name2, fid)
      temp <- paste(listName, str1, name, str2, name2, ' \n', sep='')
      cat(temp, file=fid)
      next
    }

    # For a formula
    if ("formula" %in% class(temp)) {
      temp <- paste(as.character(temp), collapse="", sep="")
      temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
      cat(temp, file=fid)
      next
    }
    # For NULL
    if (is.null(temp)) {
      temp <- paste(listName, str1, name, str2, 'NULL \n', sep='')
      cat(temp, file=fid)
      next
    }

    if (length(temp) == 1) {
      if (is.character(temp)) { 
        if (cmm == "FUNCTION") {
          temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
        } else {
          temp <- paste(listName, str1, name, str2, '"', temp, '" \n', sep='')
        }
      } else {
        temp <- paste(listName, str1, name, str2, temp, ' \n', sep='')
      }
      cat(temp, file=fid) 
    } else {
      if (flag) {
        field <- paste(listName, "$", "'", name, "'", sep="")
      } else {
        field <- paste(listName, "$", name, sep="")
      }
      #genfile.vec(temp, field, fid)
      writeToFile(temp, field, fid)
    }  
  }

  NULL

} # END: genfile.list

# Function to define a vector in swarm generator files
genfile.vec <- function(vec, name, fid) {

  if (is.null(vec)) {
    temp <- paste(name, " <- NULL\n", sep="")
    cat(temp, file=fid)
    return(NULL)
  }

  n       <- length(vec)
  vecType <- "character"
  cflag   <- 1
  if (is.numeric(vec)) {
    vecType <- "numeric"
    cflag   <- 0
  }
  temp <- paste(name, " <- ", vecType, "(", n, ") \n", sep="")
  cat(temp, file=fid)
  k <- 1
  for (temp in vec) { 
    if (cflag) {
      temp <- paste(name, '[', k, '] <- "', temp, '" \n', sep='')
    } else {
      temp <- paste(name, '[', k, '] <- ', temp, ' \n', sep='')
    }
    cat(temp, file=fid)
    k <- k + 1
  }
  
  NULL

} # END: genfile.vec


my_LocalDecode <- function(data, init, tran, pi_0, pi_1, random_draw, 
                           freq_vec, r_0, p_00, p_11, epsilon) {

  # Use same args as in EM call to reuse C code already written
  tmp       <- get_args(data, freq_vec, r_0, p_00, p_11, pi_0, pi_1, epsilon)
  iargs     <- tmp$iargs
  dargs     <- tmp$dargs
  kvec      <- numOfStates(data)
  klen      <- length(kvec)
  nr        <- nrow(data)
  nc        <- ncol(data)
  ret_vec   <- rep(-1, nr*nc)

  MISSING   <- 9
  datavec   <- as.vector(t(data))
  tmp       <- !is.finite(datavec)
  if (any(tmp)) datavec[tmp] <- MISSING

  tmp <- .C("MY_C_LocalDecode", as.integer(iargs), as.numeric(dargs), 
            as.integer(datavec), as.integer(freq_vec), as.integer(kvec),
            as.integer(random_draw), as.numeric(init), as.numeric(tran),
            ret_vec=as.integer(ret_vec))
  ret <- matrix(data=tmp$ret_vec, nrow=nr, ncol=nc, byrow=TRUE)

  ret

}

myDecode <- function(initial_data, state_pattern, init, tran, pi_0, pi_1,random_draw,
                     freq_vec, r_0, p_00, p_11, epsilon) {
  
  #Random draw decoding doesnt use patterns since we want it to be probabilistic
  #Don't need to go back to original data size since we never use the patterned version for this
  if (random_draw){
    decoded_data <- my_LocalDecode(initial_data,init,tran,pi_0,pi_1,random_draw,
                                   freq_vec, r_0, p_00, p_11, epsilon)
  } else {
    decoded_data <- myDecode0(initial_data, state_pattern, init, tran, pi_0, pi_1, random_draw,
                              freq_vec, r_0, p_00, p_11, epsilon)
  }
  return(decoded_data)
}

myPasteCols <- function(x, sep="*") {

  nc  <- ncol(x)
  ret <- x[, 1, drop=TRUE]
  if (nc > 1) {
    for (i in 2:nc) ret <- paste(ret, x[, i, drop=TRUE], sep=sep)
  }
  ret
 
}

myDecode0 <- function(initial_data, state_pattern, init, tran, pi_0, pi_1, random_draw,
                      freq_vec, r_0, p_00, p_11, epsilon){

  #Decodes data using output parameters from EM
  #locally_decoded <- LocalDecode(state_pattern,init,tran,pi_0,pi_1, random_draw)
  locally_decoded  <- my_LocalDecode(state_pattern,init,tran,pi_0,pi_1,random_draw,
                                     freq_vec, r_0, p_00, p_11, epsilon)

  #Goes back to original data size (no pattern)
  ind <- match(myPasteCols(initial_data), myPasteCols(state_pattern))
  ret <- locally_decoded[ind, , drop=FALSE]
  
  ret
}

combineSimFiles <- function(dirs, out.dir, prefix="sim22") {

  len   <- nchar(out.dir)
  if (substr(out.dir, len-1, len) != "/") out.dir <- paste0(out.dir, "/")
 
  strs  <- c("DUMY1", "DUMY2", "DUMY3", "DUMY5_boot")
  nstrs <- length(strs)
  for (i in 1:nstrs) {
    str  <- strs[i]
    out  <- paste0(out.dir, prefix, "_", str)
    fid  <- NULL
    flag <- 0
    if (i < 4) {
      ndirs <- 1
    } else {
      ndirs <- length(dirs)
    }

    for (j in 1:ndirs) {
print(c(str, j))
      ff  <- paste0(dirs[j], prefix, "_job", 1:999999, "_", str)
      tmp <- file.exists(ff)
      ff  <- ff[tmp]
      nf  <- length(ff)
      if (nf) {
        if (is.null(fid)) fid <- file(out, "w")
        for (k in 1:nf) {
          x <- scan(ff[k], sep="\n", what="character", quiet=TRUE)
          if (flag) x <- x[-1]
          if (length(x)) {
            write(x, file=fid, ncolumns=1)
            flag <- 1
          } 
        } 
      }
    }
    if (!is.null(fid)) close(fid)
  }

  NULL
}

# One dir for each model
combineSimFiles2 <- function(dirs, out.str, nRepPerFile=1) {
 
  dirs <- paste0(dirs, "/", sep="")

  ndirs <- length(dirs)
  strs  <- c("DUMY1", "DUMY2", "DUMY3", "DUMY4_boot","DUMY5_boot")
  nstrs <- length(strs)
  nrep1 <- nRepPerFile + 1

  for (i in 1:nstrs) {
    str  <- strs[i]
    out  <- paste0(out.str, "_", str)
    fid  <- file(out, "w")
    flag <- 0
    for (j in 1:ndirs) {
      dir <- dirs[j]
print(dir)
      ff  <- paste0(dir, "job", 1:99999, "_", str)
      tmp <- file.exists(ff)
      ff  <- ff[tmp]
      nf  <- length(ff)
      
      for (k in 1:nf) {
        x  <- scan(ff[k], sep="\n", what="character", quiet=TRUE)
        nx <- length(x)
        if (nRepPerFile && (nx > nrep1)) x <- x[1:nrep1]
        if (flag) x <- x[-1]
        if (length(x)) write(x, file=fid, ncolumns=1)
        flag <- 1
      } 
    }
    close(fid)
  }

  NULL
}


summarizeDUMY1 <- function(in.file, out.file) {

  cols          <- c(1:2, 4:7)
  x             <- read.table(in.file, header=1, sep="\t", stringsAsFactors=FALSE)
  nc            <- ncol(x)
  ids           <- x[, 1]
  uids          <- unique(ids)
  nids          <- length(uids)
  cx            <- colnames(x)
  bv            <- cx[-cols]
  tmp           <- !grepl("NA.", bv, fixed=TRUE)
  bvars         <- bv[tmp]
  svars         <- paste0("SD.", bvars)
  nvars         <- paste0("N.", bvars)
  nv            <- length(bvars)
  tmp           <- c(cx[cols], bvars, svars)
  ret           <- matrix(data=NA, ncol=length(tmp), nrow=nids)
  tmp[1]        <- gsub("rep_boot_", "", tmp[1], fixed=TRUE)
  colnames(ret) <- tmp

  for (v in bvars) x[, v] <- as.numeric(x[, v])
  for (i in 1:nids) {
    tmp          <- ids == uids[i]
    x2           <- x[tmp, , drop=FALSE]
    ret[i, cols] <- unlist(x2[1, cols, drop=TRUE])
    for (j in 1:nv) {
      vec              <- x2[, bvars[j], drop=TRUE]
      #ret[i, nvars[j]] <- sum(is.finite(vec)) 
      ret[i, bvars[j]] <- mean(vec, na.rm=TRUE)
      ret[i, svars[j]] <- sd(vec, na.rm=TRUE)
    }
  }
  ret[, 1] <- gsub("rep1_boot0_", "", ret[, 1], fixed=TRUE)
  if (length(out.file)) write.table(ret, file=out.file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  ret
}

summarizeDUMY2 <- function(in.file, out.file) {

  cols          <- c(1:2, 4:7)
  x             <- read.table(in.file, header=1, sep="\t", stringsAsFactors=FALSE)
  nc            <- ncol(x)
  ids           <- x[, 1]
  uids          <- unique(ids)
  nids          <- length(uids)
  cx            <- colnames(x)
  bv            <- cx[-cols]
  tmp           <- !grepl("NA.", bv, fixed=TRUE)
  bvars         <- bv[tmp]
  tmp           <- grepl("covp.", bvars, fixed=TRUE)
  covpvars      <- bvars[tmp]
  bvars         <- bvars[!tmp]
  svars         <- paste0("SD.", bvars)
  nvars         <- paste0("N.", bvars)
  nv            <- length(bvars)
  tmp           <- c(cx[cols], bvars, covpvars, svars)
  ret           <- matrix(data=NA, ncol=length(tmp), nrow=nids)
  tmp[1]        <- gsub("rep_boot_", "", tmp[1], fixed=TRUE)
  colnames(ret) <- tmp

  for (v in c(bvars, covpvars)) x[, v] <- as.numeric(x[, v])
  for (i in 1:nids) {
    tmp          <- ids == uids[i]
    x2           <- x[tmp, , drop=FALSE]
    ret[i, cols] <- unlist(x2[1, cols, drop=TRUE])
    for (j in 1:nv) {
      vec                 <- x2[, bvars[j], drop=TRUE]
      ret[i, bvars[j]]    <- mean(vec, na.rm=TRUE)
      ret[i, svars[j]]    <- sd(vec, na.rm=TRUE)

      vec                 <- x2[, covpvars[j], drop=TRUE]
      N                   <- sum(is.finite(vec)) 
      ret[i, covpvars[j]] <- sum(vec, na.rm=TRUE)/N
    }
  }
  ret[, 1] <- gsub("rep1_boot0_", "", ret[, 1], fixed=TRUE)
  if (length(out.file)) write.table(ret, file=out.file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  ret
}

summarizeDUMY3 <- function(in.file, out.file) {

  cols          <- 1:4
  x             <- read.table(in.file, header=1, sep="\t", stringsAsFactors=FALSE)
  nc            <- ncol(x)
  ids           <- x[, 1]
  uids          <- unique(ids)
  nids          <- length(uids)
  cx            <- colnames(x)
  bv            <- cx[-cols]
  tmp           <- !grepl("NA.", bv, fixed=TRUE)
  bvars         <- bv[tmp]
  svars         <- paste0("SD.", bvars)
  nvars         <- paste0("N.", bvars)
  nv            <- length(bvars)
  tmp           <- c(cx[cols], bvars, svars)
  ret           <- matrix(data=NA, ncol=length(tmp), nrow=nids)
  tmp[1]        <- gsub("rep_boot_", "", tmp[1], fixed=TRUE)
  colnames(ret) <- tmp

  for (v in bvars) x[, v] <- as.numeric(x[, v])
  for (i in 1:nids) {
    tmp          <- ids == uids[i]
    x2           <- x[tmp, , drop=FALSE]
    ret[i, cols] <- unlist(x2[1, cols, drop=TRUE])
    for (j in 1:nv) {
      vec              <- x2[, bvars[j], drop=TRUE]
      #ret[i, nvars[j]] <- sum(is.finite(vec)) 
      ret[i, bvars[j]] <- mean(vec, na.rm=TRUE)
      ret[i, svars[j]] <- sd(vec, na.rm=TRUE)
    }
  }
  ret[, 1] <- gsub("rep1_boot0_", "", ret[, 1], fixed=TRUE)
  if (length(out.file)) write.table(ret, file=out.file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  ret
}

summarizeDUMY5 <- function(in.file, out.file) {

  cols          <- 1:6
  x             <- read.table(in.file, header=1, sep="\t", stringsAsFactors=FALSE)
  nc            <- ncol(x)
  ids           <- x[, 1]
  uids          <- unique(ids)
  nids          <- length(uids)
  cx            <- colnames(x)
  bv            <- cx[-cols]
  tmp           <- !grepl("NA.", bv, fixed=TRUE)
  bvars         <- bv[tmp]
  tmp           <- grepl("covp", bvars, fixed=TRUE)
  covpvars      <- bvars[tmp]
  bvars         <- bvars[!tmp]
  tmp           <- grepl("nboot", bvars, fixed=TRUE)
  nbootvars     <- bvars[tmp]
  svars         <- paste0("SD.", bvars)
  nvars         <- paste0("N.", bvars)
  nv            <- length(bvars)
  tmp           <- c(cx[cols], bvars, covpvars, svars)
  ret           <- matrix(data=NA, ncol=length(tmp), nrow=nids)
  tmp[1]        <- gsub("rep_", "", tmp[1], fixed=TRUE)
  colnames(ret) <- tmp

  for (v in c(bvars, covpvars)) x[, v] <- as.numeric(x[, v])
  for (i in 1:nids) {
    tmp          <- ids == uids[i]
    x2           <- x[tmp, , drop=FALSE]
    ret[i, cols] <- unlist(x2[1, cols, drop=TRUE])
    for (j in 1:nv) {
      vec                 <- x2[, bvars[j], drop=TRUE]
      ret[i, bvars[j]]    <- mean(vec, na.rm=TRUE)
      ret[i, svars[j]]    <- sd(vec, na.rm=TRUE)
      vec0                <- vec

      vec                 <- x2[, covpvars[j], drop=TRUE]
      N                   <- sum(is.finite(vec)) 
      ret[i, covpvars[j]] <- sum(vec, na.rm=TRUE)/N

      if (bvars[j] %in% nbootvars) print(summary(vec0))
    }
  }
  ret[, 1] <- gsub("rep1_", "", ret[, 1], fixed=TRUE)
  ret <- as.data.frame(ret, stringsAsFactors=FALSE)
  for (i in 2:ncol(ret)) ret[, i] <- as.numeric(ret[, i])
  if (length(out.file)) myWriteTable(ret, out.file, nm="summary") 
  ret
}

myWriteTable <- function(x, out, nm="data") {
  
  len <- nchar(out)
  if (tolower(substr(out, len-4, len)) == ".xlsx") {
    obj       <- list()
    obj[[nm]] <- x
    OutputListToExcel(out, obj)
  } else {
    write.table(x, file=out, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
  }
  NULL
}

getErrorInfo <- function(dir, prefix="sim22", job=NULL) {

  if (length(job)) {
    f0   <- paste0(dir, prefix, "_job", job, "_DUMY1_boot")
    x    <- read.table(f0, header=1, sep="\t", as.is=TRUE)
    x    <- matrix(unlist(strsplit(x[, 1, drop=TRUE], "_", fixed=TRUE)), nrow=nrow(x), ncol=6, byrow=TRUE)
    x    <- x[, 2]
    x    <- gsub("boot", "", x, fixed=TRUE)
    x    <- as.numeric(x)
    tab  <- table(x)
    bids <- names(tab) 
    tmp  <- tab < 12
    bids <- as.numeric(bids[tmp])
    print(bids)
    stop()
  }

  f0   <- paste0(dir, prefix, "_job", 1:99999, "_DUMY1_boot")
  tmp  <- file.exists(f0)
  f0   <- f0[tmp]
  n    <- length(f0)
  nerr <- 0
  v    <- "beta.ib.1"
  N0   <- 13
  N1   <- (N0 - 1)*250 + 1
  TOT  <- 3001
  minn <- TOT
  minf <- ""
  for (i in 1:n) {
    #x    <- read.table(f0[i], header=1, sep="\t", as.is=TRUE)
    #tmp  <- !is.finite(as.numeric(x[, v, drop=TRUE])) 
    #m    <- sum(tmp)
    #nerr <- nerr + m
    x     <- scan(f0[i], what="char", sep="\n", quiet=TRUE)
    nx    <- length(x)
    if (nx < minn) {
      minn <- nx
      minf <- f0[i]
    }
    if (nx != TOT) {
      m    <- 1000 - (nx-1-3)/3
      nerr <- nerr + m
    } 
  }
  print(c(nerr, TOT*250))
  print(minn)
  print(minf)

  NULL

}

OutputListToExcel <- function(filename, obj) {

  N <- length(obj)
  if (!N) return(NULL)

  filename <- trimws(filename)
  len      <- nchar(filename)
  if (len < 5) stop("ERROR: filename is not valid")
  str <- tolower(substr(filename, len-4, len))
  if (str != ".xlsx") stop("ERROR: filename must have a .xlsx extension")
  if (file.exists(filename)) file.remove(filename)

  # For backwards compatibility
  if (is.data.frame(obj) || is.matrix(obj)) obj <- list(output=obj)
  N <- length(obj)

  nms <- trimws(names(obj))
  if (!length(nms)) nms <- paste("output ", 1:N, sep="")
  tmp <- nchar(nms) < 1
  if (any(tmp)) nms[tmp] <- paste("output ", (1:N)[tmp], sep="")

  over <- TRUE
  for (i in 1:N) {
    tmp <- obj[[i]]
    if (length(tmp) && (is.data.frame(tmp) || is.matrix(tmp))) {
      rio::export(tmp, filename, which=nms[i], overwrite=over)
      over <- FALSE
    }
  }
  if (file.exists(filename)) {
    msg <- paste0("Output saved to file: ", filename, "\n")
    cat(msg)
  } else {
    warning("No output written, check the input object.")
  }  

  filename
}

combineAllBoot <- function(dirs, prefix, out, nsim=1000, TEST=0) {

  flag <- 0
  reps <- paste0("rep", 1:nsim)
  fid  <- file(out, "w")
  for (i in 1:nsim) {
    ff  <- paste0(dirs, prefix, "_job", i, "_DUMY1_boot.gz")
    tmp <- file.exists(ff)
    ff  <- ff[tmp]
    nf  <- length(ff)
    if (!nf) next
    if (TEST) nf <- TEST
    for (j in 1:nf) {
      x <- scan(ff[j], what="char", sep="\n", quiet=TRUE)
      if (flag) x <- x[-1]
      flag <- 1
      x <- gsub("rep1", reps[i], x, fixed=TRUE)
      write(x, file=fid, ncolumns=1)
    }
  }
  close(fid)
  NULL
}

getSeqsFromList <- function(inlist) {

  ncomb <- 1
  nc    <- length(inlist)
  for (i in 1:nc) ncomb <- ncomb*length(inlist[[i]])
  nn    <- names(inlist)

  mat <- matrix(NA, nrow=ncomb, ncol=nc)
  if (length(nn)) colnames(mat) <- nn
  for (j in 1:nc) {
     vec  <- inlist[[j]]
     nvec <- length(vec)
     if (j == 1) {
       m <- ncomb/nvec
     } else {
       m <- m/nvec
     }
     mat[, j] <- rep(vec, each=m)
  }

  mat

} # END: getSeqsFromList
