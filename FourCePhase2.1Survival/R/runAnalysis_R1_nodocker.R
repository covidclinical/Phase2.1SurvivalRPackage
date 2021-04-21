
runAnalysis_R1_nodocker=function(currSiteId, dir.input, dir.output){
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

if(is.rmchild==1){
  patient_num.out=dat.survival$dat.analysis.deceased[which(dat.survival$dat.analysis.deceased$age_group%in%c("00to02", "03to05", "06to11", "12to17")==1), "patient_num"]
  dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)
  dat.cls=dat.cls[dat.cls$patient_num%in%patient_num.out!=1,]
}

patient_num.out=unique(dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date>"2020-10","patient_num"])
dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)

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
"neutrophil_count")

nm.lab.LabAll=colnames(dat.survival$dat.analysis.severedeceased)
nm.lab.LabAll=nm.lab.LabAll[nm.lab.LabAll%in%c("age_group", "age_cat", "race_bin", nm.cls, nm.dem, nm.dem.new, "days_since_admission", "patient_num", "value", nm.event)!=1]
nm.lab.Lit3=c("C_reactive_protein_CRP_Normal_Sensitivity","lymphocyte_count", "D_dimer")


is.rmout=F
is.onlyCLS=F

if(is.rmout==T){
  patient_num.out=lab_outlier_fun(dir.input)
  dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)}

if(is.onlyCLS==T){
  patient_num.out=setdiff(dat.survival$dat.calendar$patient_num, dat.cls$patient_num)
  dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)}

#summary.report=summary.report.fun(dat.survival, siteid, dir.output)

survfit.glmnet=NULL
period.train="all"
period.valid="all"


### 4. cov model
cat("4. covariate model \n")

survfit.coxnet.LabCommon.DemCls.impute=survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.LabCommon, nm.cls, siteid, dir.output, 
                                                             period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=F)

#survfit.coxnet.LabCommon.DemCls.ind=survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.LabCommon, nm.cls, siteid, dir.output, 
#                                                          period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=T)


survfit.coxnet.Lit3.DemCls.impute=survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.Lit3, nm.cls, siteid, dir.output, 
                                                        period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=F)

#survfit.coxnet.Lit3.DemCls.ind=survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=nm.lab.Lit3, nm.cls, siteid, dir.output, 
#                                                     period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=T)


### 5. single lab + ind
cat("5. single lab \n")

#survfit.coxnet.LabSingle.ind=lapply(nm.lab.LabAll, function(ll) tryCatch(survfit.coxnet.R1.fun(dat.survival, nm.event, nm.dem, nm.lab.keep=ll, nm.cls, siteid, dir.output, 
#                                                                                      period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=T, include.cls=F, include.dem=F), error=function(e) NA))

#names(survfit.coxnet.LabSingle.ind)=nm.lab.LabAll

survfit.lab.baseline.rm.event0=tryCatch(survfit.lab.baseline.R1.fun(dat.survival, nm.event, t0.all=c(1:14), rm.event.baseline=T, is.bt=T),error=function(e){print(e); NA})

###########
cat("6. binary model \n")

survfit.lab.t=survfit.lab.t.R1.fun(dat.survival, LocalPatientObservations, nm.event, dir.output, nm.lab.LabAll,t0.all=c(1:14), rm.event.baseline=F, is.bt=T)
survfit.lab.t.rm.event0=survfit.lab.t.R1.fun(dat.survival, LocalPatientObservations, nm.event, dir.output, nm.lab.LabAll, t0.all=c(1:14), rm.event.baseline=T, is.bt=T)


########### transportability
## coefficient from literature
cat("7. transportability \n")
data(betahat.port, package="FourCePhase2.1Survival")

betahat.Lit3=c(0.013,-1.984, 0.112/1000)
names(betahat.Lit3)=nm.lab.Lit3

survfit.coxnet.port.Lit3=tryCatch(survfit.glmnet.coefficient.R1.fun(dat.survival, ipw=T, nm.event, nm.lab.all=nm.lab.LabAll, betahat=betahat.Lit3, nm.cls, siteid, dir.output, 
                                                                    period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=F, is.bt=T),error=function(e){print(e); NA})

for(mymodel in ls(betahat.port)){
  for(submodel in ls(betahat.port[[mymodel]])){
    betahat.port[[mymodel]][[submodel]]$VA=colMeans(do.call(rbind, betahat.port[[mymodel]][[submodel]][names(betahat.port[[mymodel]][[submodel]])%in%paste0("VA", 1:5)]),na.rm=T)
  }
}
survfit.coxnet.port.betahat=NULL
site.beta.list=ls(betahat.port$LabCommon.DemCls$impute)
if("VA"%in%site.beta.list){site.beta.list=setdiff(site.beta.list, paste0("VA",1:5))}
site.europe.list=c("FRBDX", "ICSM", "APHP", "H120", "UKFR")
site.us.list=setdiff(site.beta.list, site.europe.list)

if(currSiteId%in%site.europe.list){site.beta.keep=site.us.list}else{
  site.beta.keep=intersect(site.europe.list, site.beta.list)
}

mymodel="Lit3.DemCls"
submodel="impute"

for(mysite in site.beta.keep){
  print(mysite)
  betahat=betahat.port[[mymodel]][[submodel]][[mysite]]
  survfit.coxnet.port.betahat[[mymodel]][[submodel]][[mysite]]=tryCatch(survfit.glmnet.coefficient.R1.fun(dat.survival, ipw=T, nm.event, nm.lab.all=nm.lab.LabAll, betahat= betahat, nm.cls, siteid, dir.output, 
                                                                                                          period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=F, is.bt=T),error=function(e){print(e); NA})
}


####### C statistics
cat("C statistics \n")

survfit.cstat.LabCommon.DemCls.impute=tryCatch(survfit.cstat.R1.fun(dat.survival, nm.event, nm.lab.LabCommon, nm.cls, siteid, dir.output, 
                                                                    period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=F),error=function(e){print(e); NA})

#survfit.cstat.LabCommon.DemCls.ind=tryCatch(survfit.cstat.R1.fun(dat.survival, nm.event, nm.lab.LabCommon, nm.cls, siteid, dir.output, 
#                                                        period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=T),error=function(e){print(e); NA})


survfit.cstat.Lit3.DemCls.impute=tryCatch(survfit.cstat.R1.fun(dat.survival, nm.event, nm.lab.Lit3, nm.cls, siteid, dir.output, 
                                                               period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=F),error=function(e){print(e); NA})

#survfit.cstat.Lit3.DemCls.ind=tryCatch(survfit.cstat.R1.fun(dat.survival, nm.event, nm.lab.Lit3, nm.cls, siteid, dir.output, 
#                                                   period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=T),error=function(e){print(e); NA})


save(survfit.coxnet.LabCommon.DemCls.impute=survfit.coxnet.LabCommon.DemCls.impute,
     #survfit.coxnet.LabCommon.DemCls.ind=survfit.coxnet.LabCommon.DemCls.ind,
     survfit.coxnet.Lit3.DemCls.impute=survfit.coxnet.Lit3.DemCls.impute,
     #survfit.coxnet.Lit3.DemCls.ind=survfit.coxnet.Lit3.DemCls.ind,
     #survfit.coxnet.LabSingle.ind=survfit.coxnet.LabSingle.ind,              
     survfit.lab.t=survfit.lab.t,
     survfit.lab.t.rm.event0=survfit.lab.t.rm.event0,
     survfit.lab.baseline.rm.event0=survfit.lab.baseline.rm.event0,
     survfit.coxnet.port.Lit3=survfit.coxnet.port.Lit3,
     survfit.coxnet.port.betahat=survfit.coxnet.port.betahat,
     survfit.cstat.LabCommon.DemCls.impute=survfit.cstat.LabCommon.DemCls.impute,
     #survfit.cstat.LabCommon.DemCls.ind=survfit.cstat.LabCommon.DemCls.ind,
     survfit.cstat.Lit3.DemCls.impute=survfit.cstat.Lit3.DemCls.impute,
     #survfit.cstat.Lit3.DemCls.ind=survfit.cstat.Lit3.DemCls.ind,
     file=file.path(dir.output, paste0(currSiteId, "_Result_R1.Rdata")))
}


