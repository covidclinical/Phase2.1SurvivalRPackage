
survfit.glmnet.R1.fun=function(dat.survival, ipw=T, nm.event, nm.lab.keep, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv, K=10){
  
  cat("1. data preparing \n")
  dat0=dat.prep.R1.fun(dat.survival, nm.event, nm.lab.keep, nm.cls)
  patient.rm=unique(dat0[which(dat0$days_since_admission==0& dat0$severedeceased==1),"patient_num"])
  dat0=dat0[dat0$patient_num%in%patient.rm!=1,]
  
  if(period.train=="all"){dat.train0=dat0}
  if(period.train=="early"){dat.train0=dat0[which(dat0$calendar_date<calendar.date.cut), ]}
  if(period.train=="late"){dat.train0=dat0[which(dat0$calendar_date>=calendar.date.cut), ]}
  
  if(period.valid=="all"){dat.valid0=dat0}
  if(period.valid=="early"){dat.valid0=dat0[which(dat0$calendar_date<calendar.date.cut), ]}
  if(period.valid=="late"){dat.valid0=dat0[which(dat0$calendar_date>=calendar.date.cut), ]}
  
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severedeceased", "deceased","calendar_date","calendar_day",        
                                     "age","sex","race"))
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("obs_",nm.lab.new)])
  
  multi.formulas = as.formula(paste('Y~',
                                    paste(nm.dem.new, collapse="+"),"+",
                                    paste(nm.lab.new,collapse="+")))
  
  dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(Y=0,dat.train0))[,-1])
  dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(Y=0,dat.valid0))[,-1])
  
  cat("2. training \n")
  junk.train=Est.Survfit.GLMNET(dat.train, nm.event, multi.formulas, t0.all, ipw=T)
  betahat=junk.train$betahat
  lamhat=junk.train$lamhat
  score.app=junk.train$score.app
  
  cat("3. validation \n")
  ##if there is no overlap between dat.train and dat.valid, then choose yes.cv=F
  junk.cv= CV.Survfit.GLMNET(dat.train, dat.valid, betahat, lamhat, t0.all, nm.event, ipw=T, K=10, yes.cv=T)
  score.cv=junk.cv$score.cv
  
  cat("4. accuracy parameters by days since admission\n")
  roc.cv=ROC.Survfit.FUN(dat.label=dat.valid, score=score.cv, t0.all, nm.event, ipw=T)
  
  return(list(betahat=betahat, roc.cv=roc.cv))
}
