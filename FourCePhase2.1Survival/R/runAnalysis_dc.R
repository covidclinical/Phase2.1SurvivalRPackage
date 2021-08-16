
runAnalysis_dc=function(currSiteId){
  dir.input=FourCePhase2.1Data::getInputDataDirectoryName()
  dir.output=getProjectOutputDirectory()
  obfuscation.level=FourCePhase2.1Data::getObfuscation(toupper(currSiteId))
  obfuscation=F
  if(length(obfuscation.level)!=0){
    if(obfuscation.level!=0){obfuscation=T}}
  #### 1. read data
  cat("1. read data\n")
  LocalPatientObservations=FourCePhase2.1Data::getLocalPatientObservations_nodocker(currSiteId, dir.input)
  LocalPatientClinicalCourse=FourCePhase2.1Data::getLocalPatientClinicalCourse_nodocker(currSiteId, dir.input)
  LocalPatientSummary=FourCePhase2.1Data::getLocalPatientSummary_nodocker(currSiteId, dir.input)
  data(code.dict, package="FourCePhase2.1Data")
  data(betahat.port, package="FourCePhase2.1Survival")
  data(lab.breaks.log, package="FourCePhase2.1Survival")
  data(lab.breaks.original, package="FourCePhase2.1Survival")
  data(autoimmune.icd, package="FourCePhase2.1Survival")
  
  dat=getDCData(dir.input, code.dict, currSiteId, LocalPatientClinicalCourse, LocalPatientObservations, LocalPatientSummary)
  nm.lab.keep=c("ALT","albumin", "AST", "AA","creatinine", "CRP",  "TB",  "WBC", "LYM", "neutrophil_count", "DD","charlson_score" )#"DD" ,
  nm.med.keep=c("MED.COAGB", "MED.DIURETIC","MED.REMDESIVIR","MED.ACEI", "MED.ARB" , "MED.COAGA", "MED.SIANES", "MED.SICARDIAC")
  
  nm.lab.keep=nm.lab.keep[nm.lab.keep%in%colnames(dat)]
  nm.med.keep=nm.med.keep[nm.med.keep%in%colnames(dat)]
  
  model.list=c("all-wo-autoimmune", "all-w-autoimmune", "autoimmune")
  res.DC=vector(mode="list",length=length(model.list))
  for(mymodel in model.list){
    print(mymodel)
    tryCatch({
      if(mymodel=="all-wo-autoimmune"){
        junk=paste0(paste0('Surv(days_since_admission,','deceased',')~'),paste0(c(colnames(dat)[6:8],nm.lab.keep,nm.med.keep), collapse="+") )
        multi.formulas = as.formula(junk)
        data= data.frame(dat[,1:5],model.matrix(multi.formulas, data.frame(dat))[,-1])
      }else{
        if(mymodel=="all-w-autoimmune"){  
          junk=paste0(paste0('Surv(days_since_admission,','deceased',')~'),paste0(c(colnames(dat)[6:8],nm.lab.keep,nm.med.keep, "ind_autoimmune"), collapse="+") )
          multi.formulas = as.formula(junk)
          data= data.frame(dat[,1:5],model.matrix(multi.formulas, data.frame(dat))[,-1])
        }else{
          junk=paste0(paste0('Surv(days_since_admission,','deceased',')~'),paste0(c(colnames(dat)[6:8],nm.lab.keep,nm.med.keep), collapse="+") )
          multi.formulas = as.formula(junk)
          data= data.frame(dat[dat$ind_autoimmune==1,1:5],model.matrix(multi.formulas, data.frame(dat[dat$ind_autoimmune==1,]))[,-1])  
        }
      }
      
      # further clean
      data$age00to25=(data$age00to25+data$age26to49)
      data$age26to49=NULL
      colnames(data)=gsub("age00to25", "age00to49", colnames(data))
      
      data$raceasian=(data$raceasian+data$racehispanic.and.other)
      data$racehispanic.and.other=NULL
      colnames(data)=gsub("raceasian", "raceasian.hispanic.and.other", colnames(data))
      
      # # which(apply(data[, grepl("obs_", colnames(data))==1],2,mean)<0.05)
      # # obs_cardiac_morm           obs_DD 
      # # 16               17 
      # data$cardiac_morm=NULL 
      # data$DD=NULL
      # data$obs_cardiac_morm=NULL
      # data$obs_DD=NULL
      
      # covariates=colnames(data[12:ncol(data)])
      # save(covariates,file="covariates.rda")
      # load("covariates.rda")
      # ind=match(nm.lab.keep,colnames(data[,12:ncol(data)]) )
      #data=data.frame(data[,1:11], data[,c(nm.lab.keep,nm.med.keep)])
      
      
      # analysis
      ind.test=which( (1:nrow(data)) %% 3 ==0)
      ind.train=which( (1:nrow(data)) %% 3 >0)
      data.test=data[ind.test,]
      data.train=data[ind.train,]
      
      sum_stat = surv_local(time=data.train$days_since_admission/100, event=data.train$deceased,
                            x=data.train[,6:ncol(data.train)], calendar_time=data.train$calendar_time,
                            num_bh_knot=2, d_phi=3)
      
      sum_stat$X <- NULL
      sum_stat$cbhhat.obj <- NULL
      sum_stat$alasso <- NULL
      
      sum_stat$covariates.names=colnames((data.train[,6:ncol(data.train)]))
      sum_stat$calendercutoff="2020-03-1"
      sum_stat$calender=data.frame("mean"=mean(data.train$calendar_time),"sd"=sd(data.train$calendar_time))
      x=data.train[,6:ncol(data.train)]
      sum_stat$xmean=apply(x,2,mean)
      sum_stat$xcov=cov(x)
      res.DC[[mymodel]]=sum_stat
    },error=function(e) print(e))
  }
  
  
  save(res.DC=res.DC,
       res.outcome.dist=res.outcome.dist,
       file=file.path(dir.output, paste0(currSiteId, "_DC.Rdata")))
  
}


