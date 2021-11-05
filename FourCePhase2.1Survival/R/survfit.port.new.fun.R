survfit.port.new.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, betahat, dir.output, t0.all, period.train="all", period.valid="all", method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, is.calendar=1){
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
  
  
  #dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(dat.train0))[,-1])
  #dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(dat.valid0))[,-1])
  dat.train=dat.train[,c(colnames(dat.train)[1:3],intersect(colnames(dat.train), names(betahat)))]
  dat.valid=dat.valid[,c(colnames(dat.valid)[1:3],intersect(colnames(dat.valid), names(betahat)))]
  
  cat("3. validation \n")
  betahat=betahat[names(betahat)%in%colnames(dat.valid)]
  ##if there is no overlap between dat.train and dat.valid, then choose yes.cv=F
  junk.cv.cov= CV.Survfit.COXNET.sep(dat.valid, dat.valid, betahat, lamhat=NULL, t0.all=c(1:14), nm.event, K=NULL, is.bt=F, yes.cv=F)
  score.cv=junk.cv.cov$score.sep.cov.cv
 
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

