
runAnalysis_TemporalTrend_nodocker=function(currSiteId, dir.input, dir.output){
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

#### 2. set up parameters
cat("2. set up parameters\n")

nm.cls="charlson_score"
nm.event="deceased"
is.rmout=F
is.onlyCLS=F

nm.lab.all=c("total_bilirubin","creatinine", "Ferritin", "D_dimer", "C_reactive_protein_CRP_Normal_Sensitivity", "albumin","neutrophil_count")
nm.lab.keep=c("D_dimer", "C_reactive_protein_CRP_Normal_Sensitivity", "albumin")

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

#### 4. summary report
cat("4. summary report\n")

summary.report=summary.report.fun(dat.survival, siteid, dir.output)

#### 5. coxnet
cat("5. coxnet \n")

survfit.coxnet=survfit.coxridge=NULL
period.valid="all"
nm.event="deceased"
period.train="all"
survfit.coxnet[[nm.event]][[period.train]][[period.valid]]=survfit.coxnet.fun(dat.survival, nm.event=nm.event, nm.lab.keep,nm.cls, siteid, dir.output, 
                                          period.train, period.valid, calendar.date.cut="2020-07",  t0.all=28, yes.cv=T, K=10)

survfit.coxridge[[nm.event]][[period.train]][[period.valid]]=survfit.coxridge.fun(dat.survival, nm.event=nm.event, nm.lab.keep,nm.cls, siteid, dir.output, 
                                                                                       period.train, period.valid, calendar.date.cut="2020-07",  t0.all=28, yes.cv=T, K=10)




#### 6. lab distribution
cat("6. lab distribution\n")

lab.dist.original=lab_dist_fun(dir.input, code.dict,LocalPatientObservations,dat.survival,calendar.date.cut="2020-07", myscale="original")
lab.dist.log=lab_dist_fun(dir.input,code.dict,LocalPatientObservations, dat.survival,calendar.date.cut="2020-07", myscale="log")

#### 7. lab observation rate
cat("7. lab observation rate \n")
lab.summary=lab_summary_va_fun(dir.input, code.dict, LocalPatientObservations, dat.survival)

#### 8. lab recover rate
cat("8. lab recovery rate \n")

myscale="log"
lab.recover=NULL
strat="max_day"
pat.days.cut1=0
pat.days.cut2=Inf
days.range.list=list(c(1:14))
for(iii in length(days.range.list)){
  days.range=days.range.list[[iii]]
  tmp=lab_recoverRate_new_strat_fun(dir.input, strat, dat.cls, pat.days.cut1, pat.days.cut2, days.range=days.range, code.dict, LocalPatientObservations, LocalPatientClinicalCourse, dat.survival,calendar.date.cut="2020-07", myscale, is.ns=1)
  lab.recover[[strat]][[as.character(pat.days.cut1)]][[as.character(pat.days.cut2)]][[paste(min(days.range), max(days.range), sep="-")]]=tmp
}


#### 9. charlson score
cat("9. charlson score \n")

cls=data.frame(dat.cls)
cls$patient_num=as.numeric(cls$patient_num)
dat.calendar=dat.survival$dat.calendar
junk=left_join(dat.calendar, cls, by="patient_num")
junk$calendar_date=substr(junk$calendar_date,1,7)
junk2=junk

tmp=setDT(junk)
tmp=data.frame(tmp[,lapply(.SD,mean,na.rm=TRUE),by=calendar_date])
tmp=tmp[order(tmp$calendar_date),c("calendar_date", "charlson_score")]
cls.summary=tmp

junk2$obs=1-is.na(junk2$charlson_score)
tmp=setDT(junk2)
tmp1=data.frame(tmp[,lapply(.SD,mean,na.rm=TRUE),by=calendar_date])
tmp1=tmp1[order(tmp1$calendar_date),c("calendar_date", "obs")]

tmp2=data.frame(tmp[,lapply(.SD,sum,na.rm=TRUE),by=calendar_date])
tmp2=tmp2[order(tmp2$calendar_date),c("calendar_date", "obs")]

cls.obs.summary=cbind(tmp1, tmp2[,2])
colnames(cls.obs.summary)[2]="obs_rate"
colnames(cls.obs.summary)[3]="obs_sum"

cls.early=hist(junk2$charlson_score[junk2$calendar_date<"2020-07"],plot=F)
cls.late=hist(junk2$charlson_score[junk2$calendar_date>="2020-07"], plot=F)

#### 10. obfuscation
cat("10. obfuscation\n")
if(obfuscation==T){
  junk=obfuscation.fun(summary.report, survfit.coxnet,lab.dist.original, lab.dist.log, lab.summary, lab.recover, cls.obs.summary,cls.early, cls.late, obfuscation.level)
  summary.report=junk$summary.report
  survfit.coxnet=junk$survfit.coxnet
  lab.dist.original=junk$lab.dist.original
  lab.dist.log=junk$lab.dist.log
  lab.summary=junk$lab.summary
  lab.recover=junk$lab.recover
}


save(summary.report,survfit.coxnet,survfit.coxridge,
     lab.dist.original, lab.dist.log,
     lab.recover, 
     lab.summary, 
     cls.summary, cls.obs.summary, cls.early, cls.late,file=file.path(dir.output, paste0(currSiteId, "_Result.Rdata")))

}


