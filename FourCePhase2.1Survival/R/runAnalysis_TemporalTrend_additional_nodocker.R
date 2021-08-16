
runAnalysis_TemporalTrend_additional_nodocker=function(currSiteId, dir.input, dir.output, is.maxmin=F){
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

#### 4. summary report
cat("4. summary report\n")

summary.report=summary.report.fun(dat.survival, siteid, dir.output)
dem.report=dem.report.fun(dat.survival, siteid, dir.output)

#KM <- tryCatch(survfit(Surv(days_since_admission, deceased) ~ 1, data = dat.survival$dat.analysis.deceased),error=function(e) NA)

#### 6. lab distribution
cat("6. lab distribution\n")

lab.dist.original=lab_dist_fun(dir.input, code.dict,LocalPatientObservations,dat.survival,calendar.date.cut="2020-07", myscale="original", lab.breaks.original=lab.breaks.original,lab.breaks.log=lab.breaks.log)
lab.dist.log=lab_dist_fun(dir.input,code.dict,LocalPatientObservations, dat.survival,calendar.date.cut="2020-07", myscale="log", lab.breaks.original=lab.breaks.original,lab.breaks.log=lab.breaks.log)

#### 7. lab observation rate
cat("7. lab observation rate \n")
lab.summary=lab_summary_va_fun(dir.input, code.dict, LocalPatientObservations, dat.survival)

#### 8. lab recover rate
cat("8. lab recovery rate \n")

myscale="log"
lab.recover=lab.recover.rmDead=NULL
strat="max_day"
pat.days.cut1=0
pat.days.cut2=Inf
days.range.list=list(c(1:14))
for(iii in length(days.range.list)){
  days.range=days.range.list[[iii]]
  tmp1=lab_recoverRate_new_strat_fun(dir.input, strat, dat.cls, pat.days.cut1, pat.days.cut2, days.range=days.range, nm.lab.keep=c(nm.10lab,"AA"), code.dict, LocalPatientObservations, LocalPatientClinicalCourse, dat.survival,calendar.date.cut="2020-07", myscale, is.ns=1)
  tmp2=lab_recoverRate_new_strat_rmDead_fun(dir.input, strat, dat.cls, pat.days.cut1, pat.days.cut2, days.range=days.range, nm.lab.keep=c(nm.10lab,"AA"), code.dict, LocalPatientObservations, LocalPatientClinicalCourse, dat.survival,calendar.date.cut="2020-07", myscale, is.ns=1)
  lab.recover[[strat]][[as.character(pat.days.cut1)]][[as.character(pat.days.cut2)]][[paste(min(days.range), max(days.range), sep="-")]]=tmp1
  lab.recover.rmDead[[strat]][[as.character(pat.days.cut1)]][[as.character(pat.days.cut2)]][[paste(min(days.range), max(days.range), sep="-")]]=tmp2
}

res.outcome.dist=outcome_dist_fun(LocalPatientClinicalCourse) 
patient_num.keep=unique(dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date<"2020-07","patient_num"])
res.outcome.dist.early=outcome_dist_fun(LocalPatientClinicalCourse[LocalPatientClinicalCourse$patient_num%in%patient_num.keep,]) 
patient_num.keep=unique(dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date>="2020-07","patient_num"])
res.outcome.dist.late=outcome_dist_fun(LocalPatientClinicalCourse[LocalPatientClinicalCourse$patient_num%in%patient_num.keep,]) 

save(summary.report,
     dem.report,
     lab.dist.original, lab.dist.log,
     lab.recover,lab.recover.rmDead, 
     res.outcome.dist,
     file=file.path(dir.output, paste0(currSiteId, "_TemporalTrend_additional.Rdata")))

}


