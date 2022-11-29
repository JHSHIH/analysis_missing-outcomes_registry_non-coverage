
generate.ns=function(id,times,status,ndf){
  
  ns.x=list()
  #uniq.status=c(0,1)
  uniq.status=sort(unique(status))
  for (i in uniq.status){
    ns.x[[i+1]]=list()
    temp=times[status==i]
    idtemp=id[status==i]
    for (j in 1:ndf){
      temp2=try(ns(temp,df=j),silent=T)
      if (inherits(temp2,"try-error")==F)ns.x[[i+1]][[j]]=data.frame(idtemp,temp2)
    }
  }
  ns.x
}

sample.coxph2.cmpr2=function(type,cox.model.t1,cox.model.t2,newdat.t1,newdat.t2,centime,decode,year0){
  #if USRT type<=0, consider both thyroid and non-thyroid cancer type for the failure type of registry-identified cancer
  #if USRT type==1, consider only thyroid cancer for the failure type of registry-identified cancer
  #if USRT type==2, consider only non-thyroid cancer for the failure type of registry-identified cancer
  
  captured<-rep(1,length(centime))
  impute.dat=matrix(NA,nrow(newdat.t1),2)
  nit<-50
  k<-1
  while (all(captured==0)==F&k<=nit){    
    index=which(captured==1)
    #if (length(index)>0) {
      
      decode.index<-matrix(decode[index,],nr=length(index))
      all.instate.index<-apply(decode.index==1,1,all)
      index.1<-index[all.instate.index==T]
      index.2<-index[all.instate.index==F]
      if (length(index.1)>0){impute.dat[index.1,1]<-centime[index.1]
      impute.dat[index.1,2]<-0
      captured[index.1]<-0
	}
      if (length(index.2)>0){
      year0.index<-year0[index.2]
      centime.index<-centime[index.2]

      if (type<=0){
        temp.t1=survfit(cox.model.t1,newdata=newdat.t1[index.2,], se.fit=FALSE)
        temp.t2=survfit(cox.model.t2,newdata=newdat.t2[index.2,], se.fit=FALSE)
        temp.t.surv<-matrix(temp.t1$surv*temp.t2$surv,nc=length(index.2),nr=length(temp.t1$time))
      }
      if (type==1){
        temp.t1=survfit(cox.model.t1,newdata=newdat.t1[index.2,], se.fit=FALSE)
        temp.t.surv<-matrix(temp.t1$surv,nc=length(index.2),nr=length(temp.t1$time))
      }  
      if (type==2){
        temp.t2<-survfit(cox.model.t2,newdata=newdat.t2[index.2,], se.fit=FALSE)
        temp.t.surv<-matrix(temp.t2$surv,nc=length(index.2),nr=length(temp.t2$time))
        
      }
      if (type<2)temp.t.time<-temp.t1$time else temp.t.time<-temp.t2$time
      rv1=runif(nrow(newdat.t1[index.2,]))
      newt=unlist(lapply(1:ncol(temp.t.surv),function(i)temp.t.time[which(temp.t.surv[,i]<=rv1[i])[1]]))
      newt[is.na(newt)]=centime.index[is.na(newt)]
      impute.t=apply(cbind(newt,centime.index),1,min,na.rm=T)
      impute.year<-ceiling(impute.t/12)
      impute.status=rep(0,length(newt))
      
      #index2<-which(newt<centime.index)
      index22<- which(newt<centime.index)
      temp.table<-table(newt<centime.index)
      if (length(index22)>0){
        if (type>0 ) impute.status[index22]<-type else {
          
          temp.chaz.t1<-matrix(temp.t1$cumhaz,nc=length(index),nr=length(temp.t1$time))
          temp.chaz.t2<-matrix(temp.t2$cumhaz,nc=length(index),nr=length(temp.t1$time))
          temp.haz.t1<-temp.chaz.t1-rbind(rep(0,ncol(temp.chaz.t1)),as.matrix(temp.chaz.t1[-nrow(temp.chaz.t1),]))
          temp.haz.t2<-temp.chaz.t2-rbind(rep(0,ncol(temp.chaz.t2)),as.matrix(temp.chaz.t2[-nrow(temp.chaz.t2),]))
          p.type<-temp.haz.t1/(temp.haz.t1+temp.haz.t2)
          temp.index<-match(index22,1:ncol(p.type))
          #cat(c(sum(is.na(temp.index)),temp.table,dim(p.type),length(index22)),append=T)
          #write.table(matrix(index22,nr=100),file="index22.txt",sep="\t")

	  p.type2<-t(matrix(p.type[,index22],nc=length(index22)))
          rt2<-newt[index22]
          index3<-match(rt2,temp.t.time)
          p.type3<-unlist(lapply(1:length(index3),function(i)p.type2[i,index3[i]]))
          
          impute.status[index22]<-ifelse(runif(length(index22))<=p.type3,1,2)
        }
      }
      captured.index<-unlist(lapply(1:length(index.2),function(i)ifelse(decode[index.2[i],-1][impute.year[i]+year0[index.2[i]]]==1&impute.status[i]!=0,1,0)))
      #cat(length(captured.index),sum(captured.index==0,na.rm=T),sum(is.na(captured.index)),sum(is.na(p.type3)),"&",append=T)
      captured.index[is.na(captured.index)]<-1
      impute.dat[index.2,1][captured.index==0]=impute.t[captured.index==0]
      impute.dat[index.2,2][captured.index==0]=impute.status[captured.index==0]
      
      captured[index.2]<-captured.index
      
      } #if (length(index.2)>0)
      k<-k+1
      #cat(sum(captured==1),k,"/",append=T)
  } #end while loop
      if (k>nit){
       # impute.dat[index.2,1][captured.index==1]<-centime.index[captured.index==1]
       # impute.dat[index.2,2][captured.index==1]<-0
         impute.dat[captured==1,1]<-centime[captured==1]
         impute.dat[captured==1,2]<-0
        
      }
        
  impute.dat
}  

find.quantile=function(dist,p){
  if(is.list(dist)==T)lapply(1:length(dist),function(i)which(dist[[i]]<=p[i])[1]) else
    which(dist<=p)[1]
}    

cox.prob.type<-function(temp.fit,impute.type.i){   
  temp.chaz<-temp.haz<-list()
  index.type<-match(1:length(temp.fit),impute.type.i)
  if (!is.na(index.type[1])){temp.chaz[[1]]<-matrix(temp.fit[[1]]$cumhaz,nc=length(index),nr=length(temp.fit[[1]]$time))
  temp.haz[[1]]<-temp.chaz[[1]]-rbind(rep(0,ncol(temp.chaz[[1]])),as.matrix(temp.chaz[[1]][-nrow(temp.chaz[[1]]),]))
  }
  if (!is.na(index.type[2])){temp.chaz[[2]]<-matrix(temp.fit[[2]]$cumhaz,nc=length(index),nr=length(temp.fit[[2]]$time))
  temp.haz[[2]]<-temp.chaz[[2]]-rbind(rep(0,ncol(temp.chaz[[2]])),as.matrix(temp.chaz[[2]][-nrow(temp.chaz[[2]]),]))
  } 
  if (!is.na(index.type[3])){temp.chaz[[3]]<-matrix(temp.fit[[3]]$cumhaz,nc=length(index),nr=length(temp.fit[[3]]$time))
  temp.haz[[3]]<-temp.chaz[[3]]-rbind(rep(0,ncol(temp.chaz[[3]])),as.matrix(temp.chaz[[3]][-nrow(temp.chaz[[3]]),]))
  } 
  temp.haz.all<-temp.haz[[impute.type.i[1]]]
  for (k in impute.type.i[-1])temp.haz.all<-temp.haz.all+temp.haz[[impute.type.i[k]]]
  for (k in impute.type.i)p.type[[k]]<-temp.haz[[impute.type.i[k]]]/temp.haz.all
  
  p.type
}


sample.coxph.V2.cmpr2=function(type,cox.model.t1,cox.model.t2,newdat.t1,newdat.t2,futime1,futime2,decode,year0){

  captured<-rep(1,length(futime1))
  impute.dat=matrix(NA,nrow(newdat.t1),2)
  if (is.vector(decode)==T) decode<-matrix(decode,nr=1)
    nit<-50
    k<-1
    while (all(captured==0)==F&k<=nit){    
      index=which(captured==1)
      #if (length(index)>0) {
      
      decode.index<-matrix(decode[index,],nr=length(index))
      all.instate.index<-apply(decode.index==1,1,all)
      index.1<-index[all.instate.index==T]
      index.2<-index[all.instate.index==F]
      if (length(index.1)>0){impute.dat[index.1,1]<-futime2[index.1]
      impute.dat[index.1,2]<-0
      captured[index.1]<-0
	}
      if (length(index.2)>0){
      year0.index<-year0[index.2]
      futime12<-futime1[index.2]
      futime22<-futime2[index.2]
      
      if (type<=0){
        temp.t1=survfit(cox.model.t1,newdata=newdat.t1[index.2,], se.fit=FALSE)
        temp.t2=survfit(cox.model.t2,newdata=newdat.t2[index.2,], se.fit=FALSE)
        temp.t.surv<-matrix(temp.t1$surv*temp.t2$surv,nc=length(index.2))
      }
      if (type==1){
        temp.t1=survfit(cox.model.t1,newdata=newdat.t1[index.2,], se.fit=FALSE)
        temp.t.surv<-matrix(temp.t1$surv,nc=length(index.2))
      }  
      if (type==2){
        temp.t2<-survfit(cox.model.t2,newdata=newdat.t2[index.2,], se.fit=FALSE)
        temp.t.surv<-matrix(temp.t2$surv,nc=length(index.2))
      }
      
      if (type<2)temp.t.time<-temp.t1$time else temp.t.time<-temp.t2$time
      surv.t<-unlist(lapply(1:length(futime12),function(i)temp.t.surv[max(which(temp.t.time<=futime12[i])),i]))
      rv1=runif(nrow(newdat.t1[index.2,]))*surv.t
      
      newt=unlist(lapply(1:ncol(temp.t.surv),function(i)temp.t.time[which(temp.t.surv[,i]<=rv1[i])[1]]))
      newt[is.na(newt)]=futime22[is.na(newt)]
      
      
      impute.t=apply(cbind(newt,futime22),1,min,na.rm=T)
      impute.year<-ceiling(impute.t/12)
      impute.status=rep(0,length(impute.t))
      
      index2<-which(newt<futime22)
      
      if (length(index2)>0){
        if (type>0)impute.status[index2]<-type else {
          temp.chaz.t1<-matrix(temp.t1$cumhaz,nr=length(temp.t1$time))
          temp.chaz.t2<-matrix(temp.t2$cumhaz,nr=length(temp.t2$time))
          temp.haz.t1<-temp.chaz.t1-rbind(rep(0,ncol(temp.chaz.t1)),as.matrix(temp.chaz.t1[-nrow(temp.chaz.t1),]))
          temp.haz.t2<-temp.chaz.t2-rbind(rep(0,ncol(temp.chaz.t2)),as.matrix(temp.chaz.t2[-nrow(temp.chaz.t2),]))
          p.type<-temp.haz.t1/(temp.haz.t1+temp.haz.t2)
          p.type2<-t(matrix(p.type[,index2],nc=length(index2)))
          rt2<-newt[index2]
          index3<-match(rt2,temp.t1$time)
          if (length(index2)>1)p.type3<-diag(p.type2[,index3]) else p.type3<-p.type2[,index3]
          
          impute.status[index2]<-ifelse(runif(length(index2))<=p.type3,1,2)
        }
      }
      captured.index<-unlist(lapply(1:length(index.2),function(i)ifelse(decode[index.2[i],-1][impute.year[i]+year0.index[i]]==1&impute.status[i]!=0,1,0)))
      captured.index[is.na(captured.index)]<-1

      impute.dat[index.2,1][captured.index==0]=impute.t[captured.index==0]
      impute.dat[index.2,2][captured.index==0]=impute.status[captured.index==0]
      captured[index.2]<-captured.index
      k<-k+1
      } #if(length(index.2)>0)  
    } 
    if (k>nit){
    impute.dat[captured==1,1]<-futime2[captured==1]
    impute.dat[captured==1,2]<-0
    }
    impute.dat
}  



impute.method.1.cmpr2=function(type,dat,fit.event,ftimes,decode,year0){
 #if (is.null(type)==T & delta.obs>0) only impute failure time and delta.true=delta.obs
 # if (is.null(type)==F|delta.obs==0), impute failure time and failure type
   
  #only impute failure time 
  #no mis-classfication of event type: if delta.obs=k, delta.true=k, k=1,2
  uniq.status=sort(unique(dat$delta.obs))
  impute.dat=matrix(NA,nrow(dat),2)
  captured<-rep(1,nrow(dat))
  nit<-50
  
  for (i in uniq.status){
    k<-1
    captured.index<-NULL
    repeat {    
      index=which(dat$delta.obs==i&captured==1)
      if (length(index)>0) {
        
        decode.index<-matrix(decode[index,],nr=length(index))
        all.instate.index<-apply(decode.index==1,1,all)
        index.1<-index[all.instate.index==T]
        index.2<-index[all.instate.index==F]
        if (length(index.1)>0){impute.dat[index.1,1]<-dat[index.1,3]
        impute.dat[index.1,2]<-0
        captured[index.1]<-0
	}
        if (length(index.2)>0){
        year0.index<-year0[index.2]
        subset=dat[index.2,1:2]
        cen.true<-dat[index.2,3]
        temp.surv=predict(fit.event[[i+1]],times=ftimes,newdata=subset,type="surv")
        #rv1=runif(length(temp.surv))
        rv1<-runif(length(index.2))
        index.rt=unlist(find.quantile(temp.surv,rv1))
        rt=ftimes[index.rt]
        temp=apply(cbind(rt,cen.true),1,min,na.rm=T)
        impute.t=temp
        impute.year<-ceiling(impute.t/12)
        impute.status=rep(0,nrow(subset))
        index2<-which(rt<cen.true&!is.na(index.rt))
        if (length(index2)>0){
            if (is.null(type)& i>0)impute.status[index2]<-subset$delta.obs[index2] else {
            temp.cinc<-predict(fit.event[[i+1]],times=ftimes,newdata=subset,type="risk")
            temp.inc<-list()
            for (j in 1:length(temp.cinc)){
              if (is.list(temp.cinc[[j]])==T)
               temp.inc[[j]]<-do.call(rbind,lapply(1:length(temp.cinc[[j]]),function(k)temp.cinc[[j]][[k]]-c(0,temp.cinc[[j]][[k]][1:(length(temp.cinc[[j]][[k]])-1)]))) else
               temp.inc[[j]]<-matrix(temp.cinc[[j]]-c(0,temp.cinc[[j]][1:(length(temp.cinc[[j]])-1)]),nr=1)
            }  
            p.type<-temp.inc[[1]]/(temp.inc[[1]]+temp.inc[[2]])
            #cat(i,nrow(subset),dim(p.type),max(index2,na.rm=T), max(index.rt,na.rm=T),sum(is.na(index2)),sum(is.na(index.rt)),"//",append=T)
            #if (dim(p.type)[1]<max(index2,na.rm=t))cat(summary(p.type),"??",append=T)
            if (dim(p.type)[1]>=max(index2,na.rm=t)){
	    if (length(index2)>1) p.type2<-diag(p.type[index2,index.rt[index2]]) else p.type2<-p.type[index2,index.rt[index2]]
            impute.status[index2]<-ifelse(runif(length(index2))<=p.type2,1,2)
            } else {impute.status[index2]<-0
            impute.t[index2]<-cen.true[index2]}  
	} 
        } # if (length(index)>0)
        captured.index<-unlist(lapply(1:length(index.2),function(i)ifelse(decode[index.2[i],-1][impute.year[i]+year0.index[i]]==1&impute.status[i]!=0,1,0)))
        impute.dat[index.2,1][captured.index==0]=impute.t[captured.index==0]#only impute if not captured
        impute.dat[index.2,2][captured.index==0]=impute.status[captured.index==0]
        captured[index.2]<-captured.index
        k<-k+1
        } #if(length(index.2)>0)
        if (all(captured.index==0)==T) break else {
          if (k>nit){
          impute.dat[index.2,1][captured.index==1]<-cen.true[captured.index==1]
          impute.dat[index.2,2][captured.index==1]<-0
          break}
        }
      } else break #length(index>0)
    }#repeat
    }  
  impute.dat
}  

prob.type<-function(fit.event.i,ftimes,subset,impute.type.i){
  temp.cinc<-predict(fit.event.i,times=ftimes,newdata=subset,type="risk")
  temp.inc<-p.type<-list()
  #for (j in 1:length(temp.cinc)){
  for (j in impute.type.i){
    if (is.list(temp.cinc[[j]])==T)
      temp.inc[[j]]<-do.call(rbind,lapply(1:length(temp.cinc[[j]]),function(k)temp.cinc[[j]][[k]]-c(0,temp.cinc[[j]][[k]][1:(length(temp.cinc[[j]][[k]])-1)]))) else
        temp.inc[[j]]<-matrix(temp.cinc[[j]]-c(0,temp.cinc[[j]][1:(length(temp.cinc[[j]])-1)]),nr=1)
  }  
  temp.inc.all<-temp.inc[[impute.type.i[1]]]
  for (j in impute.type.i[-1])temp.inc.all<-temp.inc.all+temp.inc[[j]]
  for (j in impute.type.i)p.type[[j]]<-temp.inc[[j]]/temp.inc.all
  p.type
}


impute.method.V.cmpr2=function(type,dat,futime1,futime2,fit.event,ftimes,decode,year0){
  uniq.status=sort(unique(dat$delta.obs))
  impute.dat=matrix(NA,nrow(dat),2)
  captured<-rep(1,nrow(dat))
  nit<-50
  
  for (i in uniq.status){
    k<-1
    captured.index<-NULL
    repeat {
      index=which(dat$delta.obs==i&captured==1)
      if (length(index)>0) {
        decode.index<-matrix(decode[index,],nr=length(index))
        all.instate.index<-apply(decode.index==1,1,all)
        index.1<-index[all.instate.index==T]
        index.2<-index[all.instate.index==F]
        if (length(index.1)>0){impute.dat[index.1,1]<-futime2[index.1]
        impute.dat[index.1,2]<-0
        captured[index.1]<-0
	}
        if (length(index.2)>0){
        year0.index<-year0[index.2]
        subset=dat[index.2,]
        subset.futime1=futime1[index.2]
        subset.futime2=futime2[index.2]
        surv.t.futime1<-diag(do.call(cbind,predict(fit.event[[i+1]],times=subset.futime1,newdata=subset,type="surv")))
        temp.surv=predict(fit.event[[i+1]],times=ftimes,newdata=subset,type="surv")
        
        rv1=runif(length(temp.surv))*surv.t.futime1
        index.rt=unlist(find.quantile(temp.surv,rv1))
        rt=ftimes[index.rt]
        temp=apply(cbind(rt,subset.futime2),1,min,na.rm=T)
        impute.t=temp
        impute.year<-ceiling(impute.t/12)
        impute.status=rep(0,nrow(subset))
        index2<-which(rt<subset.futime2&!is.na(rt))
        if (length(index2)>0){
          if (is.null(type)==T&i>0) impute.status[index2]<-i else {
            temp.cinc<-predict(fit.event[[i+1]],times=ftimes,newdata=subset,type="risk")
            temp.inc<-list()
            for (j in 1:length(temp.cinc)){
              if (is.list(temp.cinc[[j]])==T)
                temp.inc[[j]]<-do.call(rbind,lapply(1:length(temp.cinc[[j]]),function(k)temp.cinc[[j]][[k]]-c(0,temp.cinc[[j]][[k]][1:(length(temp.cinc[[j]][[k]])-1)]))) else
                  temp.inc[[j]]<-matrix(temp.cinc[[j]]-c(0,temp.cinc[[j]][1:(length(temp.cinc[[j]])-1)]),nr=1)
            }  
            
            p.type<-temp.inc[[1]]/(temp.inc[[1]]+temp.inc[[2]])
            if (length(index2)>1) p.type2<-diag(p.type[index2,index.rt[index2]]) else p.type2<-p.type[index2,index.rt[index2]]
            #p.type2<-matrix(p.type[index2,],nr=length(index2))
            #rt2<-rt[index2]
            #if (length(index2)>1)p.type3<-diag(p.type2[,rt2]) else p.type3<-p.type2[,rt2]
            impute.status[index2]<-ifelse(runif(length(index2))<=p.type2,1,2)
          }
        }
        captured.index<-unlist(lapply(1:length(index.2),function(i)ifelse(decode[index.2[i],-1][impute.year[i]+year0.index[i]]==1&impute.status[i]!=0,1,0)))
        
        impute.dat[index.2,1][captured.index==0]=impute.t[captured.index==0]
        impute.dat[index.2,2][captured.index==0]=impute.status[captured.index==0]
 
        captured[index.2]<-captured.index
        k<-k+1
        } #if(length(index.2)>0)
        if (all(captured.index==0)==T) break else {
          if (k>nit){
          impute.dat[index.2,1][captured.index==1]<-subset.futime2[captured.index==1]
          impute.dat[index.2,2][captured.index==1]<-0
          break}
        }  
      } #if (length(index)>0)
     
    } #repeat
  } 
  impute.dat
}  



surv.decode3=function(id.V,t.true.V,delta.true.V,year.dx.true.V,decoded.states.V,year.decode,year0){
  surv.time.V=delta.V=rep(NA,nrow(decoded.states.V))
  year.lastin<-unlist(lapply(1:nrow(decoded.states.V),function(i)which(decoded.states.V[i,(year0[i]+2):max(year.decode)]==0)[1]-1+year0[i]))
  #lastin=apply(as.matrix(decoded.states.V[,-1]),1,function(z)which(z==0)[1])-1
  #year.lastin=year.decode[lastin]
  surv.time.V<-t.true.V
  delta.V<-delta.true.V
  surv.time.V[!is.na(year.lastin)]<-apply(cbind((year.lastin-year0)[!is.na(year.lastin)]*12,t.true.V[!is.na(year.lastin)]),1,min)
  delta.V[!is.na(year.lastin)]<-ifelse((year.lastin-year0)[!is.na(year.lastin)]*12<t.true.V[!is.na(year.lastin)],0,delta.true.V[!is.na(year.lastin)])
  
  list(surv.time.V=surv.time.V,delta.V=delta.V)
}  


