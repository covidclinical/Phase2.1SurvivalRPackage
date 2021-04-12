survfit.cstat.R1.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv=T, K=10, is.bt=T, include.lab=T, include.dem=T, include.cls=T, include.ind=T){
  
  #cat("1. data preparing \n")
  nm.dem=c("age_group_new", "sex", "race_new")
  nm.dem.new=c("age", "sex", "race")
  if(include.ind==T){dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls)}
  if(include.ind==F){dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls); dat0=dat0[,grepl("obs_", colnames(dat0))!=1]}
  
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
  
  dat.sd=data.frame(dat0[,1:3],model.matrix(multi.formulas, data.frame(dat0))[,-1])
  dat.sd=apply(dat.sd[,-c(1:3)],2,sd,na.rm=T)
  dat.sd=dat.sd[grepl("ns.calendar",names(dat.sd))!=1]
 
  junk=CV.Survfit.CSTAT.sep(dat.train, dat.valid, nm.event, K, is.bt, yes.cv=T)
  cstat.cv=junk$cstat.cv
  cstat.se=sd(junk$cstat.bt,na.rm=T)
  
  return(list(cstat.cv=cstat.cv, 
              cstat.se=cstat.se,
              dat.sd=dat.sd
  ))
}

