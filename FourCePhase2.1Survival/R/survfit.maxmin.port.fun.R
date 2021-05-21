survfit.maxmin.port.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, beta.maxmin.int, dir.output, t0.all, include.ind=F, include.dem=T, include.lab=T, include.cls=T){
  nm.dem=c("age_group_new", "sex", "race_new")
  
  #cat("1. data preparing \n")
  if(include.ind==T){dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls)}
  if(include.ind==F){dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls); dat0=dat0[,grepl("obs_", colnames(dat0))!=1]}
  
  dat0=dat0[dat0$days_since_admission!=0,]
  dat.train0=dat0
 
  dat.valid0=dat0
 
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severe", "severedeceased", "deceased","calendar_date","calendar_day",        
                                      "age","sex","race", "charlson_score", "obs_charlson_score"))
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("obs_",nm.lab.new)])
  nm.dem.new=c("age", "sex", "race")
  
  if(sum(dat0$obs_charlson_score==0)==0){nm.cls.new=c("charlson_score")}else{
    nm.cls.new=c("charlson_score", "obs_charlson_score")}
  if(nm.event=="severedeceased"){nm.event.print="Severe-Free Survival"}
  if(nm.event=="deceased"){nm.event.print="Overall Survival"}
  if(nm.event=="severe"){nm.event.print="Non-severe Rate"}
  junk=paste0(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    ifelse(include.dem==T,paste0(paste(nm.dem.new, collapse="+"),"+"), ""),
    ifelse(include.lab==T,paste0(paste(nm.lab.new,collapse="+"),"+"),""),
    ifelse(include.lab==T & include.ind==T,paste0(paste(paste0("obs_", nm.lab.new),collapse="+"),"+"),""),
    ifelse(include.cls==T,paste0(paste(nm.cls.new,collapse="+")),""))
  if(substr(junk,nchar(junk), nchar(junk))=="+"){junk=substr(junk, 1, nchar(junk)-1)}
  multi.formulas = as.formula(junk)
  
  
  dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(dat.train0))[,-1])
  dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(dat.valid0))[,-1])
  dat.train=dat.train[,c(colnames(dat.train)[1:3],intersect(colnames(dat.train), rownames(beta.maxmin.int)))]
  dat.valid=dat.valid[,c(colnames(dat.valid)[1:3],intersect(colnames(dat.valid), rownames(beta.maxmin.int)))]
  
  cat("3. validation \n")
  ##if there is no overlap between dat.train and dat.valid, then choose yes.cv=F
  junk.cv.cov= CV.Survfit.COXNET.sep(dat.valid, dat.valid, beta.maxmin.int, lamhat=NULL, t0.all=c(1:14), nm.event, K=NULL, is.bt=F, yes.cv=F)
  score.cv=junk.cv.cov$score.cv
 
  cat("4. accuracy parameters by days since admission\n")
  roc.cov=ROC.Survfit.FUN(dat.label=dat.valid, score=score.cv, t0.all=c(1:14), nm.event, ipw=F, is.sep=T)
  
  cat("5. accuracy parameters by calendar month\n")
  #roc.cv.bymonth=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=F)
  roc.cov.bymonth=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.cv, t0.all=c(1:14), nm.event, ipw=F, is.combine=F, is.sep=T)
  
  cat("6. score categories by calendar month\n")
  #score.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.cv, roc.cv, is.combine=F, is.sep=F)
  score.cov.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.cv, roc.cov, is.combine=F, is.sep=T)
  
  return(list(roc.cov=roc.cov,
              roc.cov.bymonth=roc.cov.bymonth,
              score.cov.bymonth=score.cov.bymonth
       
  ))
}

