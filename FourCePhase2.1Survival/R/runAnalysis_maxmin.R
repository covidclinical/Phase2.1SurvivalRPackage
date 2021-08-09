
runAnalysis_maxmin=function(currSiteId){
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
  data(betamaxmin, package="FourCePhase2.1Survival")
  
  #### 2. set up parameters
  cat("2. set up parameters\n")
  
  nm.cls="charlson_score"
  nm.event="deceased"
  is.rmout=T
  
  nm.10lab=c("alanine_aminotransferase_ALT",     
             "albumin",
             "aspartate_aminotransferase_AST",
             "creatinine",
             "C_reactive_protein_CRP_Normal_Sensitivity",
             "total_bilirubin",
             "white_blood_cell_count_Leukocytes",
             "lymphocyte_count",
             "neutrophil_count", "D_dimer")
  
  #### 3. data pivot (prepare survival data and charlson)
  cat("3. data pivot\n")
  
  dat.survival0=getSurvivalData(dir.input, code.dict, siteid=currSiteId)
  nm.10lab=nm.10lab[nm.10lab%in%colnames(dat.survival0$dat.analysis.severe)]
  
  comorb=map_charlson_codes(LocalPatientObservations)
  index_scores <- comorb[[3]]
  dat.cls0<- data.frame(index_scores %>% dplyr::select(patient_num, charlson_score))
  dat.cls0$patient_num=as.numeric(dat.cls0$patient_num)
  
  dat.survival=dat.survival0
  dat.cls=dat.cls0
  dat.survival$dat.analysis.severe=left_join(dat.survival$dat.analysis.severe, dat.cls, by="patient_num")
  dat.survival$dat.analysis.deceased=left_join(dat.survival$dat.analysis.deceased, dat.cls, by="patient_num")
  dat.survival$dat.analysis.severedeceased=left_join(dat.survival$dat.analysis.severedeceased, dat.cls, by="patient_num")
  dat.survival$dat.calendar$calendar_date=substr(dat.survival$dat.calendar$calendar_date,1,7)
  dat.survival$dat.calendar=dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date<="2021-01",]
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date<"2020-03"]="2020-03"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-04"]="2020-03"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-06"]="2020-05"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-08"]="2020-07"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-10"]="2020-09"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date>="2020-12"]="2020-11"
  
  ### remove patients<18
  patient_num.out=dat.survival$dat.analysis.deceased[which(dat.survival$dat.analysis.deceased$age_group%in%c("00to02", "03to05", "06to11", "12to17")==1), "patient_num"]
  dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)
  dat.cls=dat.cls[dat.cls$patient_num%in%patient_num.out!=1,]
  dat.survival=dem.group.change.fun(dat.survival,nm.var="age_group_new",old.group.nm="00to25",new.group.nm="18to25")
  
  
  cat("10. transportability\n")
  survfit.port=vector(mode="list", length(out))
  for(nm.beta in ls(out)){
  print(nm.beta)
  betahat.tmp=out[[nm.beta]]
  if(nm.beta!="beta.site"){
  if(nm.beta=="beta.maximin"){
  id.site=which(toupper(currSiteId)==toupper(colnames(betahat.tmp)))
  betahat=betahat.tmp[,id.site]
  betahat[is.na(betahat)]=0
  names(betahat)=rownames(betahat.tmp)
  }else{
  betahat=betahat.tmp
  betahat[is.na(betahat)]=0
  }
  survfit.port[[nm.beta]]=survfit.port.new.fun(dat.survival, nm.event, nm.lab.keep=nm.10lab, nm.cls, betahat=betahat, dir.output, t0.all, period.train="all", period.valid="all", method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, is.calendar=1)
  }else{
  for(mysiteid in colnames(betahat.tmp)){
    betahat=betahat.tmp[,mysiteid]
    betahat[is.na(betahat)]=0
    names(betahat)=rownames(betahat.tmp)  
    survfit.port[[nm.beta]][[mysiteid]]=survfit.port.new.fun(dat.survival, nm.event, nm.lab.keep=nm.10lab, nm.cls, betahat=betahat, dir.output, t0.all, period.train="all", period.valid="all", method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, is.calendar=1)
  }
  }
  }
  
  save(survfit.port=survfit.port,
       file=file.path(dir.output, paste0(currSiteId, "_maxmin.Rdata")))
  
}


