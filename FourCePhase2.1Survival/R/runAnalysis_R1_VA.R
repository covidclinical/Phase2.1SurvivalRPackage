
runAnalysis_R1_VA=function(currSiteId,dir.output){
load("data/code.dict.rda")
siteid=currSiteId
obfuscation.level=10
obfuscation=F
if(length(obfuscation.level)!=0){
  if(obfuscation.level!=0){obfuscation=T}}
#### 1. read data
cat("1. read data\n")
channel = odbcDriverConnect('driver={SQL Server};server=vhacdwrb03;Database=ORD_Cho_201803054D')
sql_query.LocalPatientObservations = paste0('SELECT * FROM ORD_Cho_201803054D.HTP.[4CE_LocalPatientObservations]')
sql_query.LocalPatientClinicalCourse = paste0('SELECT * FROM ORD_Cho_201803054D.HTP.[4CE_LocalPatientClinicalCourse]')
sql_query.LocalPatientSummary = paste0('SELECT * FROM ORD_Cho_201803054D.HTP.[4CE_LocalPatientSummary]')

LocalPatientObservations0=sqlQuery(channel, sql_query.LocalPatientObservations)
LocalPatientClinicalCourse0=sqlQuery(channel, sql_query.LocalPatientClinicalCourse)
LocalPatientSummary0=sqlQuery(channel, sql_query.LocalPatientSummary)

LocalPatientObservations=LocalPatientObservations0[LocalPatientObservations0$siteid==siteid,]
LocalPatientClinicalCourse=LocalPatientClinicalCourse0[LocalPatientClinicalCourse0$siteid==siteid,]
LocalPatientSummary=LocalPatientSummary0[LocalPatientSummary0$siteid==siteid,]

cat("data pivot \n")
dir.phase1.input=paste0("P:\\ORD_Cho_201803054D\\Hong, Chuan\\4CE_data_03032021\\", siteid, "\\")
dat.survival0=getSurvivalData_va(dir.phase1.input, LocalPatientClinicalCourse, LocalPatientObservations, LocalPatientSummary, code.dict, siteid=siteid)

data <- LocalPatientObservations
comorb=map_charlson_codes(data)
index_scores <- comorb[[3]]
dat.cls0<- data.frame(index_scores %>% dplyr::select(patient_num, charlson_score))
dat.cls0$patient_num=as.numeric(dat.cls0$patient_num)

is.rmchild=0
#### 2. adding Charlson score and combine month
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
survfit.coxnet.LabCommon.DemCls.impute=survfit.coxnet.LabAll.DemCls.impute=survfit.coxnet.Lit3.DemCls.impute=NULL
survfit.lab.t=survfit.lab.t.rm.event0=survfit.lab.baseline.rm.event0=NULL
survfit.cstat.LabCommon.DemCls.impute=survfit.cstat.Lit3.DemCls.impute=NULL
for(nm.event in c("deceased", "severedeceased")){
  cat("event: ", nm.event, "\n")
  cat("4. covariate model \n")
  nm.lab.LabAll=colnames(dat.survival[[paste0("dat.analysis.", nm.event)]])
  nm.lab.LabAll=nm.lab.LabAll[nm.lab.LabAll%in%c("age_group", "age_cat", "race_bin", nm.cls, nm.dem, nm.dem.new, "days_since_admission", "patient_num", "value", nm.event)!=1]
  cat("common model \n")
  survfit.coxnet.LabCommon.DemCls.impute[[nm.event]]=tryCatch(survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.LabCommon, nm.cls, siteid, dir.output, 
                                                                                    period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, method.impute="mice", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, include.dem=1, include.lab=1, include.cls=1),error=function(e){print(e); NA})
  cat("three model \n")
  survfit.coxnet.Lit3.DemCls.impute[[nm.event]]=tryCatch(survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.Lit3, nm.cls, siteid, dir.output, 
                                                                               period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, method.impute="mice", myscale="original", is.ind=0, is.stand=0, mice.time=5, removeALT=1, include.dem=1, include.lab=1, include.cls=1),error=function(e){print(e); NA})
  
  ### 5. single lab + ind
  cat("5. single lab \n")
  survfit.lab.baseline.rm.event0[[nm.event]]=tryCatch(survfit.lab.baseline.R1.fun(dat.survival, nm.event, nm.lab.keep=nm.lab.LabCommon, nm.dem, nm.cls, t0.all=c(1:14), rm.event.baseline=T,  is.bt=T, method.impute="mice", myscale="original", is.ind=0, is.stand=0, mice.time=5,removeALT=1),error=function(e){print(e); NA})
  
  ###########
  cat("6. binary model \n")
  survfit.lab.t[[nm.event]]=tryCatch(survfit.lab.t.R1.fun(dat.survival, LocalPatientObservations, nm.event, dir.output, nm.lab.LabCommon,t0.all=c(1:14), rm.event.baseline=F, is.bt=T),error=function(e){print(e); NA})
  survfit.lab.t.rm.event0[[nm.event]]=tryCatch(survfit.lab.t.R1.fun(dat.survival, LocalPatientObservations, nm.event, dir.output, nm.lab.LabCommon, t0.all=c(1:14), rm.event.baseline=T, is.bt=T),error=function(e){print(e); NA})
}

########### transportability
cat("7. transportability \n")
###
survfit.coxnet.port.betahat.deceased=NULL
if(is.transport==T){
  data(betahat.port.deceased, package="FourCePhase2.1Survival")
  submodel="impute"
  for(mysite in ls(betahat.port.deceased$Lit3.DemCls$impute)){
    print(mysite)
    for(mymodel in c("Lit3.DemCls", "LabCommon.DemCls")){
      betahat=betahat.port.deceased[[mymodel]][[submodel]][[mysite]]
      survfit.coxnet.port.betahat.deceased[[mymodel]][[submodel]][[mysite]]=tryCatch(survfit.glmnet.coefficient.R1.fun(dat.survival, ipw=T, nm.event="deceased", nm.lab.all=nm.lab.LabAll, betahat= betahat, nm.cls, siteid, dir.output, 
                                                                                                                       period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=F, is.bt=T),error=function(e){print(e); NA})
    }
  }
}

####### C statistics
cat("C statistics \n")
for(nm.event in c("deceased","severedeceased")){
  survfit.cstat.LabCommon.DemCls.impute[[nm.event]]=tryCatch(survfit.cstat.R1.fun(dat.survival, nm.event, nm.lab.LabCommon, nm.cls, siteid, dir.output, 
                                                                                  period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time="mice", removeALT=1, include.lab=T, include.dem=T, include.cls=T),error=function(e){print(e); NA})
  
  survfit.cstat.Lit3.DemCls.impute[[nm.event]]=tryCatch(survfit.cstat.R1.fun(dat.survival, nm.event, nm.lab.Lit3, nm.cls, siteid, dir.output, 
                                                                             period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, method.impute="zero", myscale="original", is.ind=0, is.stand=0, mice.time="mice", removeALT=1, include.lab=T, include.dem=T, include.cls=T),error=function(e){print(e); NA})
}

save(survfit.coxnet.LabCommon.DemCls.impute=survfit.coxnet.LabCommon.DemCls.impute,
     survfit.coxnet.Lit3.DemCls.impute=survfit.coxnet.Lit3.DemCls.impute,
     survfit.lab.t=survfit.lab.t,
     survfit.lab.t.rm.event0=survfit.lab.t.rm.event0,
     survfit.lab.baseline.rm.event0=survfit.lab.baseline.rm.event0,
     survfit.coxnet.port.betahat.deceased=survfit.coxnet.port.betahat.deceased,
     survfit.cstat.LabCommon.DemCls.impute=survfit.cstat.LabCommon.DemCls.impute,
     survfit.cstat.Lit3.DemCls.impute=survfit.cstat.Lit3.DemCls.impute,
     file=file.path(dir.output, paste0(currSiteId, "_R1.Rdata")))



