
survfit.glmnet.coefficient.R1.fun=function(dat.survival, ipw=T, nm.event, nm.lab.all, betahat, nm.cls, siteid, dir.output, 
                            period.train, period.valid, calendar.date.cut,  t0.all, yes.cv, K=10, is.bt=T){
  
  names(betahat)=gsub("C_reactive_protein_CRP_Normal_Sensitivity", "CRP",names(betahat))
  names(betahat)=gsub("lymphocyte_count", "LYM", names(betahat))
  names(betahat)=gsub("D_dimer", "DD", names(betahat))
  names(betahat)=gsub("alanine_aminotransferase_ALT", "ALT", names(betahat))
  names(betahat)=gsub("aspartate_aminotransferase_AST", "AST", names(betahat))
  names(betahat)=gsub("total_bilirubin", "TB", names(betahat))
  names(betahat)=gsub("lactate_dehydrogenase_LDH", "LDH", names(betahat))
  names(betahat)=gsub("prothrombin_time_PT", "PT", names(betahat))
  names(betahat)=gsub("white_blood_cell_count_Leukocytes", "WBC", names(betahat))
  names(betahat)=gsub("age_group_new", "age",names(betahat))
  names(betahat)=gsub("race_new", "race", names(betahat))
  
  #cat("1. data preparing \n")
  nm.dem=c("age_group_new", "sex", "race_new")
  nm.dem.new=c("age", "sex", "race")
  dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.all, nm.cls)
  patient.rm=unique(dat0[which(dat0$days_since_admission==0& dat0$severedeceased==1),"patient_num"])
  dat0=dat0[dat0$patient_num%in%patient.rm!=1,]
  
  if(period.train=="all"){dat.train0=dat0}
  if(period.train=="early"){dat.train0=dat0[which(dat0$calendar_date<calendar.date.cut), ]}
  if(period.train=="late"){dat.train0=dat0[which(dat0$calendar_date>=calendar.date.cut), ]}
  
  if(period.valid=="all"){dat.valid0=dat0}
  if(period.valid=="early"){dat.valid0=dat0[which(dat0$calendar_date<calendar.date.cut), ]}
  if(period.valid=="late"){dat.valid0=dat0[which(dat0$calendar_date>=calendar.date.cut), ]}
  
  nm.lab.new=setdiff(colnames(dat0),c("patient_num", "days_since_admission", "severedeceased", "deceased","calendar_date","calendar_day",        
                                     "age","sex","race", "charlson_score"))
  nm.lab.new=setdiff(nm.lab.new, nm.lab.new[grepl("obs_",nm.lab.new)])
  nm.lab.new=nm.lab.new[nm.lab.new%in%names(betahat)]
  
  nm.keep=names(betahat)
  nm.keep[grepl("age", nm.keep)]="age"
  nm.keep[grepl("sex", nm.keep)]="sex"
  nm.keep[grepl("race", nm.keep)]="race"
  nm.keep=unique(nm.keep)
  nm.keep.new=nm.keep[nm.keep%in%colnames(dat0)]
  
  
  multi.formulas = as.formula(paste('Y~',
                                    paste(nm.keep.new, collapse="+")))
  
  dat.train= data.frame(dat.train0[,1:3],model.matrix(multi.formulas, data.frame(Y=0,dat.train0))[,-1])
  dat.valid= data.frame(dat.valid0[,1:3],model.matrix(multi.formulas, data.frame(Y=0,dat.valid0))[,-1])
  
  #cat("2. training \n")
  nm.keep.model=names(betahat)[names(betahat)%in%colnames(dat.train)]
  junk.train=Est.Survfit.coefficient.GLMNET(dat.train[,c("patient_num","days_since_admission",nm.event,nm.keep.model)], betahat=betahat[nm.keep.model],nm.event,  t0.all)
  score.app=junk.train$score.app
  
  #cat("3. accuracy parameters by days since admission\n")
  roc.cv=ROC.Survfit.FUN(dat.label=dat.valid, score=score.app, t0.all, nm.event, ipw=T)
  
  #cat("4. bootstrap \n")
  if(is.bt==T){
  xx=tryCatch(do.call(rbind,lapply(1:100, function(myseed){
    set.seed(myseed)
    id.bt=sample(1:dim(dat.valid)[1], replace=T)
    junk=tryCatch(ROC.Survfit.FUN(dat.label=dat.valid[id.bt,], score=score.app[id.bt,], t0.all, nm.event, ipw=T),error=function(e) NA)
    res=tryCatch(unlist(lapply(1:14, function(ll) junk[[ll]]$auc)),error=function(e) rep(NA,14))
    res
  })), error=function(e) NA)
  roc.se=tryCatch(apply(xx,2, function(x) sd(x, na.rm=T)),error=function(e) NA)
  }else{roc.se=NA}
  return(list(betahat=betahat, roc.cv=roc.cv, roc.se=roc.se))
}
