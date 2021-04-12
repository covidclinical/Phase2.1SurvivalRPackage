
survfit.lab.baseline.R1.fun=function(dat.survival, nm.event, t0.all, rm.event.baseline=F, is.bt=T){
  
  #cat("1. data preparing \n")
  
  nm.lab.all=setdiff(colnames(dat.survival$dat.analysis.deceased),c("patient_num", "days_since_admission", "deceased", "sex", "age_group", "race", "age_cat", "race_bin", "race_new","age_group_new", "charlson_score"))
  dat0=dat.prep.lab.R1.fun(dat.survival, nm.event, nm.lab.all)
  if(rm.event.baseline==T){dat0=dat0[dat0$days_since_admission!=0,]}
  
  roc.lab=roc.lab.se=NULL
  for(nm.lab in nm.lab.all){
  dat.tmp=dat0[,c("days_since_admission", nm.event, nm.lab)]
  dat.tmp[is.na(dat.tmp[,nm.lab]),nm.lab]=median(dat.tmp[is.na(dat.tmp[,nm.lab])!=1,nm.lab])
  #dat.tmp=dat.tmp[complete.cases(dat.tmp),]
  score.tmp=matrix(rep(dat.tmp[,nm.lab], length(t0.all)), ncol=length(t0.all))
  colnames(score.tmp)=t0.all
  roc.lab[[nm.lab]]=tryCatch(ROC.Survfit.FUN(dat.label=dat.tmp, score=score.tmp, t0.all, nm.event, ipw=F, is.sep=T),error=function(e) NA)
  xx=tryCatch(do.call(rbind,lapply(1:100, function(myseed){
    set.seed(myseed)
    id.bt=sample(1:dim(dat.tmp)[1], replace=T)
    junk=ROC.Survfit.FUN(dat.label=dat.tmp[id.bt,], score=score.tmp[id.bt,], t0.all, nm.event, ipw=F, is.sep=T)
    unlist(lapply(1:14, function(ll) junk[[ll]]$auc))
  })), error=function(e) NA)
  roc.lab.se[[nm.lab]]=tryCatch(apply(xx,2, function(x) sd(x, na.rm=T)),error=function(e) NA)
  }
  return(list(roc.lab=roc.lab, roc.lab.se=roc.lab.se))
}
