
survfit.lab.fun=function(dat.survival, nm.event, dir.output, t0.all){
  
  dat.calendar=dat.survival$dat.calendar
  dat.x.raw =LocalPatientObservations

  nm.lab.keep=setdiff(colnames(data_lab_clean(dat.x.raw, code.dict, nm.patient_num, 
                                             nm.days_since_admission, nm.value, day = 0)), "patient_num")
  res=NULL
  for(nm.lab in nm.lab.keep){
  dat0=dat.prep.outcome.R1.fun(dat.survival, nm.event)
  patient.rm=dat0[which(dat0$days_since_admission==0 & dat0$severedeceased==1),"patient_num"]
  dat0=dat0[dat0$patient_num%in%patient.rm!=1,]

  myroc=NULL
  for(day in t0.all){
  dat.lab = data_lab_clean(dat.x.raw, code.dict, nm.patient_num, 
                           nm.days_since_admission, nm.value, day = day)
  x=dat.lab[,c("patient_num", nm.lab)]
  dat.tmp=left_join(dat0, x, by="patient_num")
  dat.tmp=dat.tmp[complete.cases(dat.tmp), ]
  roc.tmp=tryCatch(ROC.Est.FUN(dat.tmp$severedeceased, dat.tmp[,nm.lab], yy0=0.5, fpr=seq(0,1,0.01)),error=function(e) NA)
  if(length(roc.tmp)!=1){
  if(roc.tmp[1]<0.5){roc.tmp=ROC.Est.FUN(1-dat.tmp$severedeceased, dat.tmp[,nm.lab], yy0=0.5, fpr=seq(0,1,0.01))}
  }else{
  roc.tmp=rep(NA, 613)
  }
  myroc=cbind(myroc, roc.tmp)
  }
  colnames(myroc)=paste0("day",t0.all)
  res[[nm.lab]]=myroc
  }
  names(res)=nm.lab.keep
  return(list(res=res))
}
