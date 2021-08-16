outcome_dist_fun=function (LocalPatientClinicalCourse) 
{
dat.surv.raw=LocalPatientClinicalCourse
dat.surv.raw$discharge=1-dat.surv.raw$in_hospital
daymax=max(dat.surv.raw$days_since_admission)

#### deceased
nm.event=c("deceased")
dat.sub = dat.surv.raw[dat.surv.raw[, "days_since_admission"] <= daymax, ]
patient_num.list = sort(unique(dat.sub[, "patient_num"]))
data.surv.clean.kern = function(dat.sub, patient_num, nm.event) {
  dat.tmp = dat.sub[dat.sub[, "patient_num"] == patient_num, ]
  dat.tmp=dat.tmp[order(dat.tmp$days_since_admission),]
  dat.event.tmp = dat.tmp[which(dat.tmp[, nm.event] ==  1), ]
  
  if (dim(dat.event.tmp)[1] != 0) {
    dat.end.tmp = dat.event.tmp[which.min(dat.event.tmp[,"days_since_admission"]), c("patient_num","days_since_admission", nm.event)]
  }else {
    dat.end.tmp = dat.tmp[which.max(dat.tmp[, "days_since_admission"]),c("patient_num", "days_since_admission", nm.event)]
  }
  dat.end.tmp
}
dat.surv = do.call(rbind, lapply(patient_num.list, function(xx) data.surv.clean.kern(dat.sub,xx, nm.event)))
KM.deceased <- tryCatch(survfit(Surv(days_since_admission, deceased) ~ 1, data = dat.surv),error=function(e) {print(e);return(NA)})

#### deceased or discharged
nm.event.list=c("deceased","discharge")
dat.sub = dat.surv.raw[dat.surv.raw[, "days_since_admission"] <= daymax, ]
patient_num.list = sort(unique(dat.sub[, "patient_num"]))
data.surv.clean.kern = function(dat.sub, patient_num, nm.event.list) {
  dat.tmp = dat.sub[dat.sub[, "patient_num"] == patient_num, ]
  dat.tmp=dat.tmp[order(dat.tmp$days_since_admission),]
  if(length(nm.event.list)==1){
    dat.event.tmp = dat.tmp[which(dat.tmp[, nm.event.list] ==  1), ]}else{
      dat.event.tmp = dat.tmp[which(dat.tmp[, nm.event.list[1]] == 1|dat.tmp[, nm.event.list[2]] ==  1), ]  
    }
  
  if (dim(dat.event.tmp)[1] != 0) {
    dat.end.tmp = dat.event.tmp[which.min(dat.event.tmp[,"days_since_admission"]), c("patient_num","days_since_admission", nm.event.list)]
  }else {
    dat.end.tmp = dat.tmp[which.max(dat.tmp[, "days_since_admission"]),c("patient_num", "days_since_admission", nm.event.list)]
  }
  dat.end.tmp
}
dat.surv = do.call(rbind, lapply(patient_num.list, function(xx) data.surv.clean.kern(dat.sub,xx, nm.event.list)))
dat.surv$event=rowSums(dat.surv[,nm.event.list])
dat.surv$event=ifelse(dat.surv$event>0,1,0)
KM.deceased.discharge <- tryCatch(survfit(Surv(days_since_admission, event) ~ 1, data = dat.surv),error=function(e) {print(e);return(NA)})

#### discharged
nm.event=c("discharge")
dat.sub = dat.surv.raw[dat.surv.raw[, "days_since_admission"] <= daymax, ]
patient_num.list = sort(unique(dat.sub[, "patient_num"]))
data.surv.clean.kern = function(dat.sub, patient_num, nm.event) {
  dat.tmp = dat.sub[dat.sub[, "patient_num"] == patient_num, ]
  dat.tmp=dat.tmp[order(dat.tmp$days_since_admission),]
    dat.discharge.tmp = dat.tmp[which(dat.tmp[, nm.event] ==  1), ]
    dat.deceased.tmp = dat.tmp[which(dat.tmp[, "deceased"] ==  1), ]
    
    discharge.date=Inf; mydischarge=1
    if(dim(dat.discharge.tmp)[1]!=0){discharge.date=dat.discharge.tmp[which.min(dat.discharge.tmp[,"days_since_admission"]), "days_since_admission"]}
    
    deceased.date=Inf
    if(dim(dat.deceased.tmp)[1]!=0){deceased.date=dat.deceased.tmp[which.min(dat.deceased.tmp[,"days_since_admission"]), "days_since_admission"]}
    if(deceased.date<discharge.date){discharge.date=Inf}
    if(discharge.date==Inf){mydischarge=0}
    dat.end.tmp=data.frame(patient_num, days_since_admission=discharge.date, discharge=mydischarge)
    dat.end.tmp
}
dat.surv = do.call(rbind, lapply(patient_num.list, function(xx) data.surv.clean.kern(dat.sub,xx, nm.event)))
KM.discharge <- tryCatch(survfit(Surv(days_since_admission, discharge) ~ 1, data = dat.surv),error=function(e) {print(e);return(NA)})
list(KM.deceased=KM.deceased, KM.deceased.discharge=KM.deceased.discharge, KM.discharge=KM.discharge)
}
