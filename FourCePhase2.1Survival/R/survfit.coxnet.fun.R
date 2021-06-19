survfit.coxnet.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv=T, K=10, is.bt=T, method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time, removeALT, is.calendar=1){
  
  cat("data preparing \n")
  nm.dem=c("age_group_new", "sex", "race_new")
  nm.dem.new=c("age", "sex", "race")
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
  
  dat0=dat0[dat0$days_since_admission!=0,]
  dat0=dat0[is.na(dat0$days_since_admission)!=1,]
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severe", "severedeceased", "deceased","calendar_date","calendar_day",        
                                      "age","sex","race", "charlson_score", "mis_charlson_score"))
  nm.lab.mis= nm.lab.new[grepl("mis_",nm.lab.new)]
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("mis_",nm.lab.new)])

    
  if(sum(dat0$mis_charlson_score==1)==0){nm.cls.new=c("charlson_score")}else{
    nm.cls.new=c("charlson_score", "mis_charlson_score")}
  if(nm.event=="severedeceased"){nm.event.print="Severe-Free Survival"}
  if(nm.event=="deceased"){nm.event.print="Overall Survival"}
  if(nm.event=="severe"){nm.event.print="Non-severe Rate"}
  
  if(is.ind==T){
    multi.formulas = as.formula(paste(
      paste0('Surv(days_since_admission,', nm.event,')~'),
      ifelse(is.calendar==1,paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))',"+"),""),
      paste(nm.dem.new, collapse="+"),"+",
      paste(nm.lab.new,collapse="+"), "+",
      paste(nm.lab.mis,collapse="+"),"+",
      paste(nm.cls.new,collapse="+")))
  }else{
    multi.formulas = as.formula(paste(
      paste0('Surv(days_since_admission,', nm.event,')~'),
      paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
      paste(nm.dem.new, collapse="+"),"+",
      paste(nm.lab.new,collapse="+"), "+",
      paste(nm.cls.new,collapse="+")))
  }
  
  dat.sd=data.frame(dat0[,1:3],model.matrix(multi.formulas, data.frame(dat0))[,-1])
  dat.sd=apply(dat.sd[,-c(1:3)],2,sd,na.rm=T)
  dat.sd=dat.sd[grepl("ns.calendar",names(dat.sd))!=1]
  
  dat0.stand=dat0
  if(is.stand==T){
  dat0.stand[,nm.lab.new]=dat0.stand[,nm.lab.new]/dat.sd[nm.lab.new]}
  
  if(period.train=="all"){dat.train0=dat0.stand}
  if(period.train=="early"){dat.train0=dat0.stand[which(dat0.stand$calendar_date<calendar.date.cut), ]}
  if(period.train=="late"){dat.train0=dat0.stand[which(dat0.stand$calendar_date>=calendar.date.cut), ]}
  
  if(period.valid=="all"){dat.valid0=dat0.stand}
  if(period.valid=="early"){dat.valid0=dat0.stand[which(dat0.stand$calendar_date<calendar.date.cut), ]}
  if(period.valid=="late"){dat.valid0=dat0.stand[which(dat0.stand$calendar_date>=calendar.date.cut), ]}
  
  dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(dat.train0))[,-1])
  dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(dat.valid0))[,-1])
  
  multi.formulas.dem = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.dem.new, collapse="+")))

  dat.train.dem= data.frame(dat.train0[,1:3],model.matrix(multi.formulas.dem, data.frame(dat.train0))[,-1])
  dat.valid.dem= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.dem, data.frame(dat.valid0))[,-1])
  
  if(is.ind==T){
  multi.formulas.lab = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0.stand$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.lab.new,collapse="+"),"+",
    paste(nm.lab.mis,collapse="+")))}else{
    multi.formulas.lab = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0.stand$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.lab.new,collapse="+")))
    }
  
  dat.train.lab= data.frame(dat.train0[,1:3],model.matrix(multi.formulas.lab, data.frame(dat.train0))[,-1])
  dat.valid.lab= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.lab, data.frame(dat.valid0))[,-1])
  
  multi.formulas.cls = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
    paste(nm.cls.new,collapse="+")))
  
  dat.train.cls= data.frame(dat.train0[,1:3],model.matrix(multi.formulas.cls, data.frame(dat.train0))[,-1])
  dat.valid.cls= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.cls, data.frame(dat.valid0))[,-1])
  
  cat("training \n")
  junk.train.cov=Est.Survfit.COXNET(dat.train, nm.event, multi.formulas, t0.all)
  betahat.cov=junk.train.cov$betahat.modBIC
  
  #beta.new=betahat.cov[c("AA1", "AA2")]
  #b1.new=beta.new[1]
  #b2.new=beta.new[2]
  #b1=b2.new
  #b2=b1.new-b2.new*coef.orth[2]
  
  
  lamhat.cov=junk.train.cov$lamhat
  junk.A.cov=tryCatch(Score.A.FUN(data=dat.train[,-1],betahat=betahat.cov,rtn='Score+A'), error=function(e) NA)
  negA.cov=tryCatch(junk.A.cov$neg_info_mat, error=function(e) NA)

  junk.train.dem=Est.Survfit.COXNET(dat.train.dem, nm.event, multi.formulas.dem, t0.all)
  betahat.dem=junk.train.dem$betahat.modBIC
  lamhat.dem=junk.train.dem$lamhat
  junk.A.dem=tryCatch(Score.A.FUN(data=dat.train.dem[,-1],betahat=betahat.dem,rtn='Score+A'), error=function(e) NA)
  negA.dem=tryCatch(junk.A.dem$neg_info_mat, error=function(e) NA)
  
  junk.train.lab=Est.Survfit.COXNET(dat.train.lab, nm.event, multi.formulas.lab, t0.all)
  betahat.lab=junk.train.lab$betahat.modBIC
  lamhat.lab=junk.train.lab$lamhat
  junk.A.lab=tryCatch(Score.A.FUN(data=dat.train.lab[,-1],betahat=betahat.lab,rtn='Score+A'), error=function(e) NA)
  negA.lab=tryCatch(junk.A.lab$neg_info_mat, error=function(e) NA)
  
  junk.train.cls=Est.Survfit.COXNET(dat.train.cls, nm.event, multi.formulas.cls, t0.all)
  betahat.cls=junk.train.cls$betahat.modBIC
  lamhat.cls=junk.train.cls$lamhat
  junk.A.cls=tryCatch(Score.A.FUN(data=dat.train.cls[,-1],betahat=betahat.cls,rtn='Score+A'), error=function(e) NA)
  negA.cls=tryCatch(junk.A.cls$neg_info_mat, error=function(e) NA)
  
  cat("validation \n")
  ##if there is no overlap between dat.train and dat.valid, then choose yes.cv=F
  junk.cv.cov= tryCatch(CV.Survfit.COXNET.sep(dat.train=dat.train, dat.valid=dat.valid, betahat=betahat.cov, lamhat=lamhat.cov, t0.all=t0.all, nm.event=nm.event, K=K, is.bt=T, yes.cv=T),error=function(e) NA)
  junk.cv.dem= tryCatch(CV.Survfit.COXNET.sep(dat.train=dat.train.dem, dat.valid=dat.valid.dem, betahat=betahat.dem, lamhat=lamhat.dem, t0.all=t0.all, nm.event=nm.event, K=K, is.bt=T, yes.cv=T),error=function(e) NA)
  junk.cv.lab= tryCatch(CV.Survfit.COXNET.sep(dat.train=dat.train.lab, dat.valid=dat.valid.lab, betahat=betahat.lab, lamhat=lamhat.lab, t0.all=t0.all, nm.event=nm.event, K=K, is.bt=T, yes.cv=T),error=function(e) NA)
  junk.cv.cls= tryCatch(CV.Survfit.COXNET.sep(dat.train=dat.train.cls, dat.valid=dat.valid.cls, betahat=betahat.cls, lamhat=lamhat.cls, t0.all=t0.all, nm.event=nm.event, K=K, is.bt=T, yes.cv=T),error=function(e) NA)
  
  score.sep.cov.cv=tryCatch(junk.cv.cov$score.sep.cov.cv,error=function(e) NA)
  score.sep.dem.cv=tryCatch(junk.cv.dem$score.sep.cov.cv,error=function(e) NA)
  score.sep.lab.cv=tryCatch(junk.cv.lab$score.sep.cov.cv,error=function(e) NA)
  score.sep.cls.cv=tryCatch(junk.cv.cls$score.sep.cov.cv,error=function(e) NA)
  
  junk.cov.bt=junk.cv.cov$junk.bt
  junk.dem.bt=junk.cv.dem$junk.bt
  junk.lab.bt=junk.cv.lab$junk.bt
  junk.cls.bt=junk.cv.cls$junk.bt
  
  betahat.cov.bt=tryCatch(lapply(1:length(junk.cov.bt), function(ll) junk.cov.bt[[ll]]$betahat.bt), error=function(e) NA)
  betahat.dem.bt=tryCatch(lapply(1:length(junk.dem.bt), function(ll) junk.dem.bt[[ll]]$betahat.bt), error=function(e) NA)
  betahat.lab.bt=tryCatch(lapply(1:length(junk.lab.bt), function(ll) junk.lab.bt[[ll]]$betahat.bt), error=function(e) NA)
  betahat.cls.bt=tryCatch(lapply(1:length(junk.cls.bt), function(ll) junk.cls.bt[[ll]]$betahat.bt), error=function(e) NA)
  
  score.sep.cov.bt=tryCatch(lapply(1:length(junk.cov.bt),function(ll) junk.cov.bt[[ll]]$score.sep.cov.bt),error=function(e) NA)
  score.sep.dem.bt=tryCatch(lapply(1:length(junk.dem.bt),function(ll) junk.dem.bt[[ll]]$score.sep.cov.bt),error=function(e) NA)
  score.sep.lab.bt=tryCatch(lapply(1:length(junk.lab.bt),function(ll) junk.lab.bt[[ll]]$score.sep.cov.bt),error=function(e) NA)
  score.sep.cls.bt=tryCatch(lapply(1:length(junk.cls.bt),function(ll) junk.cls.bt[[ll]]$score.sep.cov.bt),error=function(e) NA)
  
  #score.sep.dem.bt=tryCatch(junk.cv.dem$score.sep.cov.bt,error=function(e) NA)
  #score.sep.lab.bt=tryCatch(junk.cv.lab$score.sep.cov.bt,error=function(e) NA)
  #score.sep.cls.bt=tryCatch(junk.cv.cls$score.sep.cov.bt,error=function(e) NA)
  
  cat("accuracy parameters by days since admission\n")
  roc.sep.cov.cv=tryCatch(ROC.Survfit.FUN(dat.label=dat.valid, score=score.sep.cov.cv, t0.all, nm.event, ipw=F, is.sep=T),error=function(e) NA)
  roc.sep.dem.cv=tryCatch(ROC.Survfit.FUN(dat.label=dat.valid.dem, score=score.sep.dem.cv, t0.all, nm.event, ipw=F, is.sep=T),error=function(e) NA)
  roc.sep.lab.cv=tryCatch(ROC.Survfit.FUN(dat.label=dat.valid.lab, score=score.sep.lab.cv, t0.all, nm.event, ipw=F, is.sep=T),error=function(e) NA)
  roc.sep.cls.cv=tryCatch(ROC.Survfit.FUN(dat.label=dat.valid.cls, score=score.sep.cls.cv, t0.all, nm.event, ipw=F, is.sep=T),error=function(e) NA)

  roc.sep.cov.bt=lapply(1:length(score.sep.cov.bt), function(ll) tryCatch(ROC.Survfit.FUN.bt(dat.label=dat.valid, score=score.sep.cov.bt[[ll]], t0.all,nm.event=nm.event, ipw=F, is.sep=T),error=function(e) NA))
  roc.sep.dem.bt=lapply(1:length(score.sep.dem.bt), function(ll) tryCatch(ROC.Survfit.FUN.bt(dat.label=dat.valid, score=score.sep.dem.bt[[ll]], t0.all, nm.event=nm.event, ipw=F, is.sep=T),error=function(e) NA))
  roc.sep.lab.bt=lapply(1:length(score.sep.lab.bt), function(ll) tryCatch(ROC.Survfit.FUN.bt(dat.label=dat.valid, score=score.sep.lab.bt[[ll]], t0.all, nm.event=nm.event, ipw=F, is.sep=T),error=function(e) NA))
  roc.sep.cls.bt=lapply(1:length(score.sep.cls.bt), function(ll) tryCatch(ROC.Survfit.FUN.bt(dat.label=dat.valid, score=score.sep.cls.bt[[ll]], t0.all, nm.event=nm.event, ipw=F, is.sep=T),error=function(e) NA))
  
  cat("accuracy parameters by calendar month\n")
  roc.sep.cov.cv.bymonth=tryCatch(ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.cov.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA)
  roc.sep.dem.cv.bymonth=tryCatch(ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.dem.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA)
  roc.sep.lab.cv.bymonth=tryCatch(ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.lab.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA)
  roc.sep.cls.cv.bymonth=tryCatch(ROC.Survfit.bymonth.FUN(dat.label=dat.valid0, score=score.sep.cls.cv, t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA)

  roc.sep.cov.bt.bymonth=lapply(1:length(score.sep.cov.bt), function(ll) tryCatch(ROC.Survfit.bymonth.FUN.bt(dat.label=dat.valid0, score=score.sep.cov.bt[[ll]], t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA))
  roc.sep.dem.bt.bymonth=lapply(1:length(score.sep.dem.bt), function(ll) tryCatch(ROC.Survfit.bymonth.FUN.bt(dat.label=dat.valid0, score=score.sep.dem.bt[[ll]], t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA))
  roc.sep.lab.bt.bymonth=lapply(1:length(score.sep.lab.bt), function(ll) tryCatch(ROC.Survfit.bymonth.FUN.bt(dat.label=dat.valid0, score=score.sep.lab.bt[[ll]], t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA))
  roc.sep.cls.bt.bymonth=lapply(1:length(score.sep.cls.bt), function(ll) tryCatch(ROC.Survfit.bymonth.FUN.bt(dat.label=dat.valid0, score=score.sep.cls.bt[[ll]], t0.all, nm.event, ipw=F, is.combine=F, is.sep=T),error=function(e) NA))
  
  cat("score categories by calendar month\n")
  score.sep.cov.bymonth=tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.cov.cv, roc.sep.cov.cv, is.combine=F, is.sep=T),error=function(e) NA)
  score.sep.dem.bymonth=tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.dem.cv, roc.sep.dem.cv, is.combine=F, is.sep=T),error=function(e) NA)
  score.sep.lab.bymonth=tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.lab.cv, roc.sep.lab.cv, is.combine=F, is.sep=T),error=function(e) NA)
  score.sep.cls.bymonth=tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.cls.cv, roc.sep.cls.cv, is.combine=F, is.sep=T),error=function(e) NA)

  score.sep.cov.bymonth.bt=lapply(1:length(score.sep.cov.bt),function(ll) tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.cov.bt[[ll]], roc.sep.cov.bt[[ll]], is.combine=F, is.sep=T),error=function(e) NA))
  score.sep.dem.bymonth.bt=lapply(1:length(score.sep.dem.bt),function(ll) tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.dem.bt[[ll]], roc.sep.dem.bt[[ll]], is.combine=F, is.sep=T),error=function(e) NA))
  score.sep.lab.bymonth.bt=lapply(1:length(score.sep.lab.bt),function(ll) tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.lab.bt[[ll]], roc.sep.lab.bt[[ll]], is.combine=F, is.sep=T),error=function(e) NA))
  score.sep.cls.bymonth.bt=lapply(1:length(score.sep.cls.bt),function(ll) tryCatch(score.summary(dat.label=dat.valid0, nm.event, score.sep.cls.bt[[ll]], roc.sep.cls.bt[[ll]], is.combine=F, is.sep=T),error=function(e) NA))
  
  return(list(betahat.cov=betahat.cov, 
              betahat.dem=betahat.dem, 
              betahat.lab=betahat.lab, 
              betahat.cls=betahat.cls, 
              betahat.cov.bt=betahat.cov.bt, 
              betahat.dem.bt=betahat.dem.bt, 
              betahat.lab.bt=betahat.lab.bt, 
              betahat.cls.bt=betahat.cls.bt, 
              negA.cov=negA.cov,
              negA.dem=negA.dem,
              negA.lab=negA.lab,
              negA.cls=negA.cls,
              dat.sd=dat.sd,
              roc.sep.cov.cv=roc.sep.cov.cv,
              roc.sep.dem.cv=roc.sep.dem.cv, 
              roc.sep.lab.cv=roc.sep.lab.cv,
              roc.sep.cls.cv=roc.sep.cls.cv,
              roc.sep.cov.bt=roc.sep.cov.bt,
              roc.sep.dem.bt=roc.sep.dem.bt, 
              roc.sep.lab.bt=roc.sep.lab.bt,
              roc.sep.cls.bt=roc.sep.cls.bt,
              roc.sep.cov.cv.bymonth=roc.sep.cov.cv.bymonth,
              roc.sep.dem.cv.bymonth=roc.sep.dem.cv.bymonth, 
              roc.sep.lab.cv.bymonth=roc.sep.lab.cv.bymonth,
              roc.sep.cls.cv.bymonth=roc.sep.cls.cv.bymonth,
              roc.sep.cov.bt.bymonth=roc.sep.cov.bt.bymonth,
              roc.sep.dem.bt.bymonth=roc.sep.dem.bt.bymonth, 
              roc.sep.lab.bt.bymonth=roc.sep.lab.bt.bymonth,
              roc.sep.cls.bt.bymonth=roc.sep.cls.bt.bymonth,
              score.sep.cov.bymonth=score.sep.cov.bymonth,
              score.sep.dem.bymonth=score.sep.dem.bymonth, 
              score.sep.lab.bymonth=score.sep.lab.bymonth, 
              score.sep.cls.bymonth=score.sep.cls.bymonth,
              score.sep.cov.bymonth.bt=score.sep.cov.bymonth.bt,
              score.sep.dem.bymonth.bt=score.sep.dem.bymonth.bt, 
              score.sep.lab.bymonth.bt=score.sep.lab.bymonth.bt, 
              score.sep.cls.bymonth.bt=score.sep.cls.bymonth.bt
  ))
}

