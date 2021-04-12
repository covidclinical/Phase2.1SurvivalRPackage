survfit.coxridge.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv=T, K=50){
  
  cat("data preparing \n")
  nm.dem=c("age_group_new", "sex", "race_new")
  nm.dem.new=c("age", "sex", "race")
  dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls)
  dat0=dat0[dat0$days_since_admission!=0,]
  if(period.train=="all"){dat.train0=dat0}
  if(period.train=="early"){dat.train0=dat0[which(dat0$calendar_date<calendar.date.cut), ]}
  if(period.train=="late"){dat.train0=dat0[which(dat0$calendar_date>=calendar.date.cut), ]}
  
  if(period.valid=="all"){dat.valid0=dat0}
  if(period.valid=="early"){dat.valid0=dat0[which(dat0$calendar_date<calendar.date.cut), ]}
  if(period.valid=="late"){dat.valid0=dat0[which(dat0$calendar_date>=calendar.date.cut), ]}
  
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severe", "severedeceased", "deceased","calendar_date","calendar_day",        
                                      "age","sex","race", "charlson_score", "obs_charlson_score"))
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("obs_",nm.lab.new)])
  if(sum(dat0$obs_charlson_score==0)==0){nm.cls.new=c("charlson_score")}else{
  nm.cls.new=c("charlson_score", "obs_charlson_score")}
  if(nm.event=="severedeceased"){nm.event.print="Severe-Free Survival"}
  if(nm.event=="deceased"){nm.event.print="Overall Survival"}
  if(nm.event=="severe"){nm.event.print="Non-severe Rate"}
  
  multi.formulas = as.formula(paste(
                                    paste0('Surv(days_since_admission,', nm.event,')~'),
                                    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day),15),collapse=','),'))'),"+",
                                    paste(nm.dem.new, collapse="+"),"+",
                                    paste(nm.lab.new,collapse="+"), "+",
                                    paste(paste0("obs_", nm.lab.new),collapse="+"),"+",
                                    paste(nm.cls.new,collapse="+")))
  
  dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(dat.train0))[,-1])
  dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(dat.valid0))[,-1])
  
  dat.sd=data.frame(dat0[,1:3],model.matrix(multi.formulas, data.frame(dat0))[,-1])
  dat.sd=apply(dat.sd[,-c(1:3)],2,sd,na.rm=T)
  dat.sd=dat.sd[grepl("ns.calendar",names(dat.sd))!=1]
  multi.formulas.dem = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day),15),collapse=','),'))'),"+",
    paste(nm.dem.new, collapse="+")))

  dat.train.dem= data.frame(dat.train0[,1:3],model.matrix(multi.formulas.dem, data.frame(dat.train0))[,-1])
  dat.valid.dem= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.dem, data.frame(dat.valid0))[,-1])
  
  multi.formulas.lab = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day),15),collapse=','),'))'),"+",
    paste(nm.lab.new,collapse="+"),"+",
    paste(paste0("obs_", nm.lab.new),collapse="+")))
  
  dat.train.lab= data.frame(dat.train0[,1:3],model.matrix(multi.formulas.lab, data.frame(dat.train0))[,-1])
  dat.valid.lab= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.lab, data.frame(dat.valid0))[,-1])
  
  multi.formulas.cls = as.formula(paste(
    paste0('Surv(days_since_admission,', nm.event,')~'),
    paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day),15),collapse=','),'))'),"+",
    paste(nm.cls.new,collapse="+")))
  
  dat.train.cls= data.frame(dat.train0[,1:3],model.matrix(multi.formulas.cls, data.frame(dat.train0))[,-1])
  dat.valid.cls= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas.cls, data.frame(dat.valid0))[,-1])
  
  cat("training \n")
  junk.train.cov=Est.Survfit.COXRIDGE(dat.train, nm.event, multi.formulas, t0.all)
  betahat.cov=junk.train.cov$betahat
  junk.A.cov=tryCatch(Score.A.FUN(data=dat.train[,-1],betahat=betahat.cov,rtn='Score+A'), error=function(e) NA)
  negA.cov=tryCatch(junk.A.cov$neg_info_mat, error=function(e) NA)
  
  junk.train.dem=Est.Survfit.COXRIDGE(dat.train.dem, nm.event, multi.formulas.dem, t0.all)
  betahat.dem=junk.train.dem$betahat
  junk.A.dem=tryCatch(Score.A.FUN(data=dat.train.dem[,-1],betahat=betahat.dem,rtn='Score+A'), error=function(e) NA)
  negA.dem=tryCatch(junk.A.dem$neg_info_mat, error=function(e) NA)
  
  junk.train.lab=Est.Survfit.COXRIDGE(dat.train.lab, nm.event, multi.formulas.lab, t0.all)
  betahat.lab=junk.train.lab$betahat
  junk.A.lab=tryCatch(Score.A.FUN(data=dat.train.lab[,-1],betahat=betahat.lab,rtn='Score+A'), error=function(e) NA)
  negA.lab=tryCatch(junk.A.lab$neg_info_mat, error=function(e) NA)
  
  junk.train.cls=Est.Survfit.COXRIDGE(dat.train.cls, nm.event, multi.formulas.cls, t0.all)
  betahat.cls=junk.train.cls$betahat
  junk.A.cls=tryCatch(Score.A.FUN(data=dat.train.cls[,-1],betahat=betahat.cls,rtn='Score+A'), error=function(e) NA)
  negA.cls=tryCatch(junk.A.cls$neg_info_mat, error=function(e) NA)
  
  return(list(betahat.cov=betahat.cov, 
              betahat.dem=betahat.dem, 
              betahat.lab=betahat.lab, 
              betahat.cls=betahat.cls,
              negA.cov=negA.cov,
              negA.dem=negA.dem,
              negA.lab=negA.lab,
              negA.cls=negA.cls
  ))
}

