
runAnalysis_SensitivityAnalysis=function(currSiteId){
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
data(beta.maxmin, package="FourCePhase2.1Survival")

#### 2. set up parameters
cat("2. set up parameters\n")

nm.cls="charlson_score"
nm.event="deceased"
is.rmout=F
is.onlyCLS=F

nm.lab.all=c("total_bilirubin","creatinine", "Ferritin", "D_dimer", "C_reactive_protein_CRP_Normal_Sensitivity", "albumin","neutrophil_count")
nm.3lab=c("D_dimer", "C_reactive_protein_CRP_Normal_Sensitivity", "albumin")
nm.9lab=c("alanine_aminotransferase_ALT",     
                   "albumin",
                   "aspartate_aminotransferase_AST",
                   "creatinine",
                   "C_reactive_protein_CRP_Normal_Sensitivity",
                   "total_bilirubin",
                   "white_blood_cell_count_Leukocytes",
                   "lymphocyte_count",
                   "neutrophil_count")

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

patient_num.out=dat.survival$dat.analysis.deceased[which(dat.survival$dat.analysis.deceased$age_group%in%c("00to02", "03to05", "06to11", "12to17")==1), "patient_num"]
dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)
dat.cls=dat.cls[dat.cls$patient_num%in%patient_num.out!=1,]

if(is.rmout==T){
  patient_num.out=lab_outlier_fun(dir.input)
  dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)}

if(is.onlyCLS==T){
  patient_num.out=setdiff(dat.survival$dat.calendar$patient_num, dat.cls$patient_num)
  dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)}


#### 5. coxnet
cat("5. coxnet \n")

survfit.coxnet.mice=survfit.coxnet.mice.ind=NULL
period.valid="all"
nm.event="deceased"
period.train="all"
model.setting="10lab"
print(model.setting)
nm.lab.keep=get(paste0("nm.", model.setting))
survfit.coxnet.mice[[nm.event]][[period.train]][[period.valid]][[model.setting]]=survfit.coxnet.new.fun(dat.survival, nm.event=nm.event, nm.lab.keep,nm.cls, siteid, dir.output, 
                                          period.train, period.valid, calendar.date.cut="2020-07",  t0.all=28, yes.cv=T, K=10, method.impute="mice", is.ind=F)
survfit.coxnet.mice.ind[[nm.event]][[period.train]][[period.valid]][[model.setting]]=survfit.coxnet.new.fun(dat.survival, nm.event=nm.event, nm.lab.keep,nm.cls, siteid, dir.output, 
                                                                                                    period.train, period.valid, calendar.date.cut="2020-07",  t0.all=28, yes.cv=T, K=10, method.impute="MICE", is.ind=T)
save(survfit.coxnet.mice, survfit.coxnet.mice.ind, file=file.path(dir.output, paste0(currSiteId, "_SensitivityAnalysis_Result.Rdata")))
}


