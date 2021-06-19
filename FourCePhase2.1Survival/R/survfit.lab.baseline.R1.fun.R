
survfit.lab.baseline.R1.fun=function(dat.survival, nm.event, nm.lab.keep, nm.dem, nm.cls, t0.all, rm.event.baseline=F, is.bt=T, method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time,removeALT){
  
  #cat("1. data preparing \n")
  dat00=dat.prep.mice.scale.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls, method.impute, myscale, is.ind, mice.time)
  dat0=dat00
  if("ALT"%in%colnames(dat0) & "AST"%in%colnames(dat0)){
    X.ALT=exp(dat0$ALT)
    X.AST=exp(dat0$AST)
    if(is.ind==1){
      mis_AA=ifelse(dat0$mis_ALT+dat0$mis_AST>0,1,0)
      dat0=data.frame(dat0, AA=X.AST/X.ALT, mis_AA=mis_AA)}else{
        dat0=data.frame(dat0, AA=X.AST/X.ALT)  
      }
    if(removeALT==1){
      dat0=dat0[,grepl("ALT", colnames(dat0))!=1]
    }
  }
  
  if(rm.event.baseline==T){dat0=dat0[dat0$days_since_admission!=0,]}
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severe", "severedeceased", "deceased","calendar_date","calendar_day",        
                                      "age","sex","race", "charlson_score", "mis_charlson_score"))
  nm.lab.mis= nm.lab.new[grepl("mis_",nm.lab.new)]
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("mis_",nm.lab.new)])
  
  nm.lab.all=nm.lab.new
  roc.lab=roc.lab.bt=NULL
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
  roc.lab.bt[[nm.lab]]=xx
  }
  return(list(roc.lab=roc.lab, roc.lab.bt=roc.lab.bt))
}
