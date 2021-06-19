survfit.cstat.R1.fun=function(dat.survival, nm.event, nm.lab.keep, nm.dem, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv=T, K=10, is.bt=T, method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time="mice", removeALT=1, include.lab=T, include.dem=T, include.cls=T){
  
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
  
  junk=CV.Survfit.CSTAT.sep(dat.train, dat.valid, nm.event, K, is.bt, yes.cv=T)
  cstat.cv=junk$cstat.cv
  cstat.bt=junk$cstat.bt
  
  return(list(cstat.cv=cstat.cv, 
              cstat.bt=cstat.bt,
              dat.sd=dat.sd
  ))
}

