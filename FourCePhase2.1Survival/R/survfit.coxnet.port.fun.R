survfit.coxnet.port.fun=function(dat.survival, 
                                 beta.cov, beta.dem, beta.lab, beta.cls,
                                 roc.sep.cov.cv, roc.sep.dem.cv, roc.sep.lab.cv, roc.sep.cls.cv,
                                 nm.event, nm.lab.keep, nm.cls, dir.output, t0.all){
  
  cat("1. data preparing \n")
  
  dat0=dat.prep.fun(dat.survival, nm.event, nm.lab.keep, nm.cls)
  dat0=dat0[dat0$days_since_admission!=0,]
  dat.valid0=dat0
  
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severedeceased", "deceased","calendar_date","calendar_day",        
                                      "age","sex","race", "charlson_score", "obs_charlson_score"))
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("obs_",nm.lab.new)])
  nm.dem.new=c("age", "sex", "race")
  if(sum(dat0$obs_charlson_score==0)==0){nm.cls.new=c("charlson_score")}else{
  nm.cls.new=c("charlson_score", "obs_charlson_score")}
  if(nm.event=="severedeceased"){nm.event.print="Severe-Free Survival"}
  if(nm.event=="deceased"){nm.event.print="Overall Survival"}
  if(nm.event=="severe"){nm.event.print="Non-severe Rate"}
  
  multi.formulas = as.formula(paste(
                                    paste0('Surv(days_since_admission,', nm.event,')~'),
                                    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
                                    paste(nm.dem.new, collapse="+"),"+",
                                    paste(nm.lab.new,collapse="+"), "+",
                                    paste(paste0("obs_", nm.lab.new),collapse="+"),"+",
                                    paste(nm.cls.new,collapse="+")))
  
  dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(dat.valid0))[,-1])
  
  
  multi.formulas.dem = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.dem.new, collapse="+")))

  dat.valid.dem= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.dem, data.frame(dat.valid0))[,-1])
  
  multi.formulas.lab = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.lab.new,collapse="+"),"+",
    paste(paste0("obs_", nm.lab.new),collapse="+")))
  
  dat.valid.lab= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.lab, data.frame(dat.valid0))[,-1])
  
  multi.formulas.cls = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.cls.new,collapse="+")))
  
  dat.valid.cls= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.cls, data.frame(dat.valid0))[,-1])
  
  cat("3. validation \n")
  ##if there is no overlap between dat.train and dat.valid, then choose yes.cv=F
  junk.cv.cov= CV.Survfit.COXNET.sep(dat.valid, dat.valid, betahat.cov, lamhat=NULL, t0.all, nm.event, K=NULL, yes.cv=F)
  junk.cv.dem= CV.Survfit.COXNET.sep(dat.valid.dem, dat.valid.dem, betahat.dem, lamhat=NULL, t0.all, nm.event, K=NULL, yes.cv=F)
  junk.cv.lab= CV.Survfit.COXNET.sep(dat.valid.lab, dat.valid.lab, betahat.lab, lamhat=NULL, t0.all, nm.event, K=NULL, yes.cv=F)
  junk.cv.cls= CV.Survfit.COXNET.sep(dat.valid.cls, dat.valid.cls, betahat.cls, lamhat=NULL, t0.all, nm.event, K=NULL, yes.cv=F)

  score.sep.cov.cv=junk.cv.cov$score.sep.cov.cv
  score.sep.dem.cv=junk.cv.dem$score.sep.cov.cv
  score.sep.lab.cv=junk.cv.lab$score.sep.cov.cv
  score.sep.cls.cv=junk.cv.cls$score.sep.cov.cv

  cat("4. accuracy parameters by days since admission\n")
  roc.sep.cov.cv.new=ROC.Survfit.FUN(dat.label=dat.valid, score=score.sep.cov.cv, t0.all, nm.event, ipw=F, is.sep=T)
  roc.sep.dem.cv.new=ROC.Survfit.FUN(dat.label=dat.valid.dem, score=score.sep.dem.cv, t0.all, nm.event, ipw=F, is.sep=T)
  roc.sep.lab.cv.new=ROC.Survfit.FUN(dat.label=dat.valid.lab, score=score.sep.lab.cv, t0.all, nm.event, ipw=F, is.sep=T)
  roc.sep.cls.cv.new=ROC.Survfit.FUN(dat.label=dat.valid.cls, score=score.sep.cls.cv, t0.all, nm.event, ipw=F, is.sep=T)

  cat("5. accuracy parameters by calendar month\n")
  #roc.cv.bymonth=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=F)
  roc.sep.cov.cv.bymonth.new=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.cov.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T)
  roc.sep.dem.cv.bymonth.new=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.dem.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T)
  roc.sep.lab.cv.bymonth.new=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.lab.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T)
  roc.sep.cls.cv.bymonth.new=ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.cls.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T)

  cat("6. score categories by calendar month\n")
  #score.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.cv, roc.cv, is.combine=F, is.sep=F)
  score.sep.cov.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.sep.cov.cv, roc.sep.cov.cv, is.combine=F, is.sep=T)
  score.sep.dem.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.sep.dem.cv, roc.sep.dem.cv, is.combine=F, is.sep=T)
  score.sep.lab.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.sep.lab.cv, roc.sep.lab.cv, is.combine=F, is.sep=T)
  score.sep.cls.bymonth=score.summary(dat.label=dat.valid0, nm.event, score.sep.cls.cv, roc.sep.cls.cv, is.combine=F, is.sep=T)

  return(list(roc.sep.cov.cv.new=roc.sep.cov.cv.new,
              roc.sep.dem.cv.new=roc.sep.dem.cv.new, 
              roc.sep.lab.cv.new=roc.sep.lab.cv.new,
              roc.sep.cls.cv.new=roc.sep.cls.cv.new,
              roc.sep.cov.cv.bymonth.new=roc.sep.cov.cv.bymonth.new,
              roc.sep.dem.cv.bymonth.new=roc.sep.dem.cv.bymonth.new, 
              roc.sep.lab.cv.bymonth.new=roc.sep.lab.cv.bymonth.new,
              roc.sep.cls.cv.bymonth.new=roc.sep.cls.cv.bymonth.new,
              score.sep.cov.cv=score.sep.cov.cv, 
              score.sep.dem.cv=score.sep.dem.cv, 
              score.sep.lab.cv=score.sep.lab.cv,
              score.sep.cls.cv=score.sep.cls.cv,
              score.sep.cov.bymonth=score.sep.cov.bymonth,
              score.sep.dem.bymonth=score.sep.dem.bymonth, 
              score.sep.lab.bymonth=score.sep.lab.bymonth, 
              score.sep.cls.bymonth=score.sep.cls.bymonth
  ))
}

