
runAnalysis_R1_nodocker_strat=function(currSiteId, dir.input, dir.output, is.transport=T){
obfuscation.level=getObfuscation(toupper(currSiteId))
obfuscation=F
if(length(obfuscation.level)!=0){
  if(obfuscation.level!=0){obfuscation=T}}
#### 1. read data
cat("1. read data\n")
LocalPatientObservations=FourCePhase2.1Data::getLocalPatientObservations_nodocker(currSiteId, dir.input)
LocalPatientClinicalCourse=FourCePhase2.1Data::getLocalPatientClinicalCourse_nodocker(currSiteId, dir.input)
LocalPatientSummary=FourCePhase2.1Data::getLocalPatientSummary_nodocker(currSiteId, dir.input)
data(code.dict, package="FourCePhase2.1Data")
cat("data pivot \n")

dat.survival0=getSurvivalData(dir.input, code.dict, siteid=currSiteId)
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

patient_num.out=unique(dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date>"2020-10","patient_num"])
dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)

### remove patients<18
patient_num.out=dat.survival$dat.analysis.deceased[which(dat.survival$dat.analysis.deceased$age_group%in%c("00to02", "03to05", "06to11", "12to17")==1), "patient_num"]
dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)
dat.cls=dat.cls[dat.cls$patient_num%in%patient_num.out!=1,]
dat.survival=dem.group.change.fun(dat.survival,nm.var="age_group_new",old.group.nm="00to25",new.group.nm="18to25")

#### 3. setting up variable names
nm.dem=c("age_group_new", "sex", "race_new")
nm.dem.new=c("age", "sex", "race")
nm.event="severedeceased"
nm.cls="charlson_score"

nm.lab.LabCommon=c("alanine_aminotransferase_ALT",     
           "albumin",
           "aspartate_aminotransferase_AST",
           "creatinine",
           "C_reactive_protein_CRP_Normal_Sensitivity",
           "total_bilirubin",
           "white_blood_cell_count_Leukocytes",
           "lymphocyte_count",
           "neutrophil_count", "D_dimer")

nm.lab.Lit3=c("C_reactive_protein_CRP_Normal_Sensitivity","lymphocyte_count", "D_dimer")


#summary.report=summary.report.fun(dat.survival, siteid, dir.output)

survfit.glmnet=NULL
period.train="all"
period.valid="all"

### 4. cov model
survfit.coxnet.LabCommon.DemCls.impute=survfit.coxnet.LabCommon.DemCls.impute.bt=NULL
for(nm.event in c("deceased")){
  cat("event: ", nm.event, "\n")
  cat("4. covariate model \n")
  cat("common model \n")
  survfit.coxnet.LabCommon.DemCls.impute[[nm.event]]=tryCatch(survfit.coxnet.R1.strat.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.LabCommon, nm.cls, siteid, dir.output,                                                                        period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, method.impute="mice", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, include.dem=1, include.lab=1, include.cls=1),error=function(e){print(e); NA})
}

survfit.coxnet.LabCommon.DemCls.impute.bt[[nm.event]]=tryCatch(survfit.coxnet.R1.strat.bt.fun(dat.survival.bt, nm.event, nm.dem, nm.lab.keep=nm.lab.LabCommon, nm.cls, siteid, dir.output,                                                                        period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, method.impute="mice", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, include.dem=1, include.lab=1, include.cls=1),error=function(e){print(e); NA})

save(survfit.coxnet.LabCommon.DemCls.impute=survfit.coxnet.LabCommon.DemCls.impute,
     survfit.coxnet.LabCommon.DemCls.impute.bt=survfit.coxnet.LabCommon.DemCls.impute.bt,
     file=file.path(dir.output, paste0(currSiteId, "_R1_strat.Rdata")))
}

