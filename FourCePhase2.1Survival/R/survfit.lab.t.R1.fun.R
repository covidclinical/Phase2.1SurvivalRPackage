
survfit.lab.t.R1.fun=function(dat.survival, LocalPatientObservations, nm.event, dir.output, nm.lab.keep, t0.all, rm.event.baseline=F, is.bt=T){
  
  dat.calendar=dat.survival$dat.calendar
  dat.x.raw =LocalPatientObservations

  res=NULL
  for(nm.lab in nm.lab.keep){
  print(nm.lab)
  dat0=dat.prep.outcome.R1.fun(dat.survival, nm.event)
  if(rm.event.baseline==T){
  patient.rm=dat0[which(dat0$days_since_admission==0 & dat0$severedeceased==1),"patient_num"]
  dat0=dat0[dat0$patient_num%in%patient.rm!=1,]
  }

  myroc=myroc.se=NULL
  for(day in t0.all){
  dat.lab = data_lab_clean3(dat.x.raw, code.dict, nm.value="value", day = day)
  if(nm.lab%in%colnames(dat.lab)){
  x=dat.lab[,c("patient_num", nm.lab)]
  dat.tmp=left_join(dat0, x, by="patient_num")
  dat.tmp=dat.tmp[complete.cases(dat.tmp), ]
  if(length(unique(dat.tmp[,nm.event],na.rm=T))!=1){
  roc.tmp=tryCatch(ROC.Est.FUN(dat.tmp$severedeceased, dat.tmp[,nm.lab], yy0=0.5, fpr=seq(0,1,0.01)),error=function(e) NA)
  if(is.na(roc.tmp[1])!=1){
  if(roc.tmp[1]<0.5){roc.tmp=ROC.Est.FUN(1-dat.tmp$severedeceased, dat.tmp[,nm.lab], yy0=0.5, fpr=seq(0,1,0.01))}
  }else{roc.tmp=rep(NA, 613)}
  }else{roc.tmp=rep(NA, 613)}
  }else{roc.tmp=rep(NA,613)}
  
  if(is.bt==T){roc.tmp.bt=lapply(1:100, function(myseed){
    set.seed(myseed)
    dat.bt=dat.tmp[sample(1:dim(dat.tmp)[1], dim(dat.tmp)[1],replace=T), ]
      res=tryCatch(ROC.Est.FUN(dat.bt$severedeceased, dat.bt[,nm.lab], yy0=0.5, fpr=seq(0,1,0.01))[1],error=function(e) NA)
      res
  })}else{roc.tmp.bt=NA}

  myroc=cbind(myroc, roc.tmp)
  myroc.se=c(myroc.se,sd(unlist(roc.tmp.bt),na.rm=T))
  }
  colnames(myroc)=paste0("day",t0.all)
  res[[nm.lab]]=list(myroc=myroc, myroc.se=myroc.se)
  }
  names(res)=nm.lab.keep
  return(list(res=res))
}
