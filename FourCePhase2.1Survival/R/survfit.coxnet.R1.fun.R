survfit.coxnet.R1.fun=function(dat.survival, nm.event, nm.lab.keep, nm.dem, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv=T, K=10, is.bt=T, method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time, removeALT, include.dem=1, include.lab=1, include.cls=1){
  
  #cat("1. data preparing \n")
  nm.dem=c("age_group_new", "sex", "race_new")
  nm.dem.new=c("age", "sex", "race")
  dat00=dat.prep.mice.scale.fun(dat.surv
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
    junk=paste(
      paste0('Surv(days_since_admission,', nm.event,')~'),
      ifelse(include.dem==T,paste0(paste(nm.dem.new, collapse="+"),"+"), ""),
      ifelse(include.lab==T,paste0(paste(nm.lab.new,collapse="+"),"+"),""),
      ifelse(include.lab==T,paste0(paste(nm.lab.mis,collapse="+"),"+"),""),
      ifelse(include.cls==T,paste0(paste(nm.cls.new,collapse="+")),""))
  }else{
      junk=paste(
        paste0('Surv(days_since_admission,', nm.event,')~'),
        ifelse(include.dem==T,paste0(paste(nm.dem.new, collapse="+"),"+"), ""),
        ifelse(include.lab==T,paste0(paste(nm.lab.new,collapse="+"),"+"),""),
        ifelse(include.cls==T,paste0(paste(nm.cls.new,collapse="+")),""))
  }
  if(substr(junk,nchar(junk), nchar(junk))=="+"){junk=substr(junk, 1, nchar(junk)-1)}
  multi.formulas = as.formula(junk)
  
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
  
  
  #cat("2. training \n")
  junk.train.cov=Est.Survfit.COXNET(dat.train, nm.event, multi.formulas, t0.all)
  betahat.cov=junk.train.cov$betahat.modBIC
  lamhat.cov=junk.train.cov$lamhat
  junk.A.cov=tryCatch(Score.A.FUN(data=dat.train[,-1],betahat=betahat.cov,rtn='Score+A'), error=function(e) NA)
  negA.cov=tryCatch(junk.A.cov$neg_info_mat, error=function(e) NA)

  
  #cat("3. validation \n")
  ##if there is no overlap between dat.train and dat.valid, then choose yes.cv=F
  junk.cv.cov= tryCatch(CV.Survfit.COXNET.sep(dat.train=dat.train, dat.valid=dat.valid, betahat=betahat.cov, lamhat=lamhat.cov, t0.all=t0.all, nm.event=nm.event, K=K, is.bt=T, yes.cv=T),error=function(e) NA)
  
  #cat("4. accuracy parameters by days since admission\n")
  score.sep.cov.cv=tryCatch(junk.cv.cov$score.sep.cov.cv,error=function(e) NA)
  roc.sep.cov.cv=tryCatch(ROC.Survfit.FUN(dat.label=dat.valid, score=score.sep.cov.cv, t0.all, nm.event, ipw=F, is.sep=T),error=function(e) NA)
  
  junk.cov.bt=junk.cv.cov$junk.bt
  score.sep.cov.bt=tryCatch(lapply(1:length(junk.cov.bt),function(ll) junk.cov.bt[[ll]]$score.sep.cov.bt),error=function(e) NA)
  roc.sep.cov.bt=lapply(1:length(score.sep.cov.bt), function(ll) tryCatch(ROC.Survfit.FUN.bt(dat.label=dat.valid, score=score.sep.cov.bt[[ll]], t0.all,nm.event=nm.event, ipw=F, is.sep=T),error=function(e) NA))
  
  betahat.cov.bt=tryCatch(lapply(1:length(junk.cov.bt), function(ll) junk.cov.bt[[ll]]$betahat.bt), error=function(e) NA)
  
  xx=t(data.matrix(dat.train[,-(1:3)]))%*%data.matrix(dat.train[,-(1:3)])
  N.train=dim(dat.train)[1]
  N.event=sum(dat.train[,nm.event])
  return(list(betahat.cov=betahat.cov, 
              negA.cov=negA.cov,
              dat.sd=dat.sd,
              roc.sep.cov.cv=roc.sep.cov.cv,
              roc.sep.cov.bt=roc.sep.cov.bt,
              xx=xx,
              N.train=N.train,
              N.event=N.event
  ))
}

