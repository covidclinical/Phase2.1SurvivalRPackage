
runAnalysis_dc_nodocker=function(currSiteId, dir.input, dir.output){
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
  data(betamaxmin, package="FourCePhase2.1Survival")
  data(autoimmune.icd, package="FourCePhase2.1Survival")
  
  
  dc_dat=getDCData(code.dict, LocalPatientClinicalCourse, LocalPatientObservations, LocalPatientSummary)
  dat = dc_dat$impute
  dat_nonimpute = dc_dat$non_impute
  
  nm.lab.keep=c("ALT","albumin", "AST", "AA","creatinine", "CRP",  "TB",  "WBC", "LYM", "neutrophil_count", "DD","charlson_score" )#"DD" ,
  nm.med.keep=c("MED.COAGB", "MED.DIURETIC","MED.REMDESIVIR","MED.ACEI", "MED.ARB" , "MED.COAGA")
  
  nm.lab.keep=nm.lab.keep[nm.lab.keep%in%colnames(dat)]
  nm.med.keep=nm.med.keep[nm.med.keep%in%colnames(dat)]
  
  if(sum(!is.na(dat$race))==0|length(unique(dat$race))==1){
    model.list=c("all-wo-race")
  }else{
    model.list=c("all-wo-race", "all-w-race", "white",  "black", "other", "nonwhite")
  }
  res.DC=c()
  for(mymodel in model.list){
    print(mymodel)
    tryCatch({
      data = splitdata_race(mymodel, dat, nm.lab.keep, nm.med.keep)
      data_nonimpute = splitdata_race(mymodel, dat_nonimpute, nm.lab.keep, nm.med.keep)
      
      # Calculate missing rate
      miss_rate = apply(data_nonimpute[,-c(1:5)], 2, function(x){sum(is.na(x))})/nrow(data)
      drop_var = names(miss_rate[miss_rate > 0.98])
      if(length(drop_var)!=0){
        data = data[ ,!colnames(data) %in% drop_var] 
      }
      
      # further clean
      if(!is.null(data$age00to25) & !is.null(data$age26to49)){
        data$age00to25=(data$age00to25+data$age26to49)
        data$age26to49=NULL
        colnames(data)=gsub("age00to25", "age00to49", colnames(data))
      }else{
        if(is.null(data$age00to25) & !is.null(data$age26to49)){
          data$age00to49=data$age26to49
        }else{
          if(!is.null(data$age00to25) & is.null(data$age26to49)){
            data$age00to49=data$age00to25
          }
        }
      }

      n.cv = 3
      # analysis
      
      for(kk in 1:n.cv){
        ind.test=which( (1:nrow(data)) %% n.cv == kk-1)
        ind.train=setdiff(1:nrow(data), ind.test)
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
        res.DC[[mymodel]][[kk]]=sum_stat
      }
    },error=function(e) print(e))
  }
  dat_race_summary = table(dat$race)
  save(res.DC=res.DC,
       race_table = dat_race_summary,
       file=file.path(dir.output, paste0(currSiteId, "_DC.Rdata")))
}


