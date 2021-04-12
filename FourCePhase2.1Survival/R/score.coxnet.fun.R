score.coxnet.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all){
  
  cat("1. data preparing \n")
  
  dat0=dat.prep.fun(dat.survival, nm.event, nm.lab.keep, nm.cls)
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
  nm.dem.new=c("age", "sex", "race")
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

  cat("2. training \n")
  junk.train.cov=Est.Survfit.COXNET(dat.train, nm.event, multi.formulas, t0.all)
  betahat.cov=junk.train.cov$betahat.modBIC
  lamhat.cov=junk.train.cov$lamhat
  
  junk=do.call(cbind,lapply(t0.all, function(t0) 
    1-FUN.predict.cox(newZ=dat.valid[,-c(1:3)],beta=betahat.cov,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
  ))
  score=data.frame(patient_num=dat.valid$patient_num, score_death=junk)
  return(score)
}

