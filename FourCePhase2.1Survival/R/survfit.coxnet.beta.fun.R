survfit.coxnet.beta.fun=function(dat.survival, nm.event, nm.lab.keep, nm.cls, siteid, dir.output, myscale, is.ind=T, my.ind="obs", is.stand=F){
  
  cat("data preparing \n")
  nm.dem=c("age_group_new", "sex", "race_new")
  nm.dem.new=c("age", "sex", "race")
  dat0=dat.prep.mice.scale.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls, method.impute="zero", myscale=myscale,is.ind=T)
  dat0=dat0[dat0$days_since_admission!=0,]
  #dat0.obs=dat0[,grepl("obs_", colnames(dat0))]
  #dat0.mis=1-dat0.obs
  #colnames(dat0.mis)=gsub("obs_", "mis_", colnames(dat0.mis))
  #dat0=cbind(dat0, dat0.mis)
  dat.train0=dat0

  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severe", "severedeceased", "deceased","calendar_date","calendar_day",        
                                      "age","sex","race", "charlson_score", "mis_charlson_score"))
  nm.lab.new=setdiff(nm.lab.new, c(nm.lab.new[grepl("obs_",nm.lab.new)], nm.lab.new[grepl("mis_",nm.lab.new)]))
  if(sum(dat0$mis_charlson_score==0)==0){nm.cls.new=c("charlson_score")}else{
  nm.cls.new=c("charlson_score", "mis_charlson_score")}
  if(nm.event=="severedeceased"){nm.event.print="Severe-Free Survival"}
  if(nm.event=="deceased"){nm.event.print="Overall Survival"}
  if(nm.event=="severe"){nm.event.print="Non-severe Rate"}
  
  if(is.ind==T){
      multi.formulas = as.formula(paste(
      paste0('Surv(days_since_admission,', nm.event,')~'),
      paste0('ns(calendar_day, knots=c(',paste(seq(15,max(dat0$calendar_day)-1,15),collapse=','),'))'),"+",
      paste(nm.dem.new, collapse="+"),"+",
      paste(nm.lab.new,collapse="+"), "+",
      paste(paste0("mis_", nm.lab.new),collapse="+"),"+",
      paste(nm.cls.new,collapse="+")))
  }else{
    multi.formulas = as.formula(paste(
      paste0('Surv(days_since_admission,', nm.event,')~'),
      paste(nm.lab.new,collapse="+"),"+",
      paste(nm.dem.new, collapse="+"),"+",
      paste(nm.cls.new,collapse="+")))
  }
  
  dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(dat.train0))[,-1])
  
  dat.sd=data.frame(dat0[,1:3],model.matrix(multi.formulas, data.frame(dat0))[,-1])
  dat.sd=apply(dat.sd[,-c(1:3)],2,sd,na.rm=T)
  dat.sd=dat.sd[grepl("ns.calendar",names(dat.sd))!=1]

  dat.train.stand=dat.train
  dat.train.stand[,nm.lab.new]=dat.train.stand[,nm.lab.new]/dat.sd[nm.lab.new]
    
  if(is.stand==T){dat.train.final=dat.train.stand}else{
    dat.train.final=dat.train
  }
  cat("training \n")
  junk=Est.ALASSO.GLMNET(dat.train.final[,-1],fam0="Cox", Wi=NULL)
  mybeta=junk$bhat.modBIC
  mybeta=mybeta[which(grepl("ns.calendar_day",names(mybeta))!=1)]
  return(list(mybeta=mybeta, dat.sd=dat.sd, dat=dat.train.final))
}

