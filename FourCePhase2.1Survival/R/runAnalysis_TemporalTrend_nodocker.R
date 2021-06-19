
runAnalysis_TemporalTrend_nodocker=function(currSiteId, dir.input, dir.output, is.maxmin=F){
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

KM <- tryCatch(survfit(Surv(days_since_admission, deceased) ~ 1, data = dat.survival$dat.analysis.deceased),error=function(e) NA)

#### 5. coxnet
cat("5. coxnet \n")

survfit.coxnet=survfit.coxridge=NULL
period.valid="all"
nm.event="deceased"
period.train="all"
model.setting ="10lab"
myscale="original"
print(model.setting)
nm.lab.keep=get(paste0("nm.", model.setting))
for(removeALT in c(1)){
for(is.stand in c(0)){
for(method.impute in c("zero", "mice")){
  print(method.impute)
  if(method.impute=="mice"){is.ind.list=c(0,1); mice.time.list=c(5)}else{is.ind.list=1; mice.time.list=NA}
      for(is.ind in is.ind.list){
        print(is.ind)
        for(mice.time in mice.time.list){
        survfit.coxnet[[method.impute]][[myscale]][[paste0("ind",is.ind)]][[paste0("stand",is.stand)]][[paste0("mice", mice.time)]][[paste0("removeALT",removeALT)]]=survfit.coxnet.fun(dat.survival, nm.event=nm.event, nm.lab.keep,nm.cls, siteid, dir.output, 
                                                                                                                          period.train, period.valid, calendar.date.cut="2020-07",  t0.all=28, yes.cv=T, K=10, is.bt=T, method.impute=method.impute, myscale=myscale, is.ind=is.ind, is.stand=is.stand, mice.time=mice.time, removeALT=removeALT)
        }}
}}
}
 
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
  tmp1=lab_recoverRate_new_strat_fun(dir.input, strat, dat.cls, pat.days.cut1, pat.days.cut2, days.range=days.range, nm.lab.keep=nm.10lab, code.dict, LocalPatientObservations, LocalPatientClinicalCourse, dat.survival,calendar.date.cut="2020-07", myscale, is.ns=1)
  tmp2=lab_recoverRate_new_strat_rmDead_fun(dir.input, strat, dat.cls, pat.days.cut1, pat.days.cut2, days.range=days.range, nm.lab.keep=nm.10lab, code.dict, LocalPatientObservations, LocalPatientClinicalCourse, dat.survival,calendar.date.cut="2020-07", myscale, is.ns=1)
  lab.recover[[strat]][[as.character(pat.days.cut1)]][[as.character(pat.days.cut2)]][[paste(min(days.range), max(days.range), sep="-")]]=tmp1
  lab.recover.rmDead[[strat]][[as.character(pat.days.cut1)]][[as.character(pat.days.cut2)]][[paste(min(days.range), max(days.range), sep="-")]]=tmp2
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
tmp.mean=data.frame(tmp[,lapply(.SD,mean,na.rm=TRUE),by=calendar_date])
tmp.mean=tmp.mean[order(tmp.mean$calendar_date),c("calendar_date", "charlson_score")]
colnames(tmp.mean)[2]="charlson_score_mean"
tmp.sd=data.frame(tmp[,lapply(.SD,sd,na.rm=TRUE),by=calendar_date])
tmp.sd=tmp.sd[order(tmp.sd$calendar_date),c("calendar_date", "charlson_score")]
colnames(tmp.sd)[2]="charlson_score_sd"

cls.summary=data.frame(tmp.mean, charlson_score_sd=tmp.sd[,2])

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

cat("10. maxmin\n")
survfit.maxmin.port=NULL
if(is.maxmin==T){
  id.maxmin=which(toupper(currSiteId)==toupper(ls(beta.maxmin)))
  beta.maxmin.int=beta.maxmin[[id.maxmin]]
  survfit.maxmin.port=survfit.maxmin.port.fun(dat.survival, nm.event, nm.9lab, nm.cls, beta.maxmin.int, dir.output, t0.all, include.ind=F, include.dem=T, include.lab=T, include.cls=T)
}

#### 11. obfuscation
cat("11. obfuscation\n")
if(obfuscation==T){
  junk=obfuscation.fun(summary.report, KM, survfit.coxnet,lab.dist.original, lab.dist.log, lab.summary, lab.recover, lab.recover.rmDead, cls.obs.summary,cls.early, cls.late, obfuscation.level)
  summary.report=junk$summary.report
  KM=junk$KM
  survfit.coxnet=junk$survfit.coxnet
  lab.dist.original=junk$lab.dist.original
  lab.dist.log=junk$lab.dist.log
  lab.summary=junk$lab.summary
  lab.recover=junk$lab.recover
  lab.recover=junk$lab.recover
  
}

save(summary.report,
     dem.report,
     KM, 
     survfit.coxnet,
     lab.dist.original, lab.dist.log,
     lab.recover,lab.recover.rmDead, 
     lab.summary, 
     cls.summary, cls.obs.summary, cls.early, cls.late,
     survfit.maxmin.port=survfit.maxmin.port,
     file=file.path(dir.output, paste0(currSiteId, "_TemporalTrend.Rdata")))

}


