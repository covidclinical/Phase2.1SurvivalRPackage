lab_summary_fun=function(dir.input, dat.survival){
dat.calendar=dat.survival$dat.calendar
dat.x.raw = read.csv(paste0(dir.input, "/LocalPatientObservations.csv"))
patient_num.severe=dat.survival$dat.analysis.severe$patient_num[which(dat.survival$dat.analysis.severe$severe==1)]
patient_num.deceased=dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1)]
patient_num.severedeceased=dat.survival$dat.analysis.severedeceased$patient_num[which(dat.survival$dat.analysis.severedeceased$severedeceased==1)]
res=NULL
for(nm.event in c("all", "deceased", "severedeceased")){
  if(nm.event=="all"){
    lab.obs=lab.obs.rate.fun(dat.x.raw, dat.calendar, nm.patient_num, nm.days_since_admission)
    lab.mean=lab.obs.mean.fun(dat.x.raw, dat.calendar, nm.patient_num, nm.days_since_admission)
    lab.mean.std=lab.obs.mean.std.fun(dat.x.raw, dat.calendar, nm.patient_num, nm.days_since_admission)
    }
  if(nm.event=="severedeceased"){
    lab.obs=lab.obs.rate.fun(dat.x.raw[dat.x.raw$patient_num%in%patient_num.severedeceased,], dat.calendar[dat.calendar$patient_num%in%patient_num.severedeceased,], nm.patient_num, nm.days_since_admission)
    lab.mean=lab.obs.mean.fun(dat.x.raw[dat.x.raw$patient_num%in%patient_num.severedeceased,], dat.calendar[dat.calendar$patient_num%in%patient_num.severedeceased,], nm.patient_num, nm.days_since_admission)
    lab.mean.std=lab.obs.mean.std.fun(dat.x.raw[dat.x.raw$patient_num%in%patient_num.severedeceased,], dat.calendar[dat.calendar$patient_num%in%patient_num.severedeceased,], nm.patient_num, nm.days_since_admission) 
  }
  if(nm.event=="deceased"){
    lab.obs=lab.obs.rate.fun(dat.x.raw[dat.x.raw$patient_num%in%patient_num.deceased,], dat.calendar[dat.calendar$patient_num%in%patient_num.deceased,], nm.patient_num, nm.days_since_admission)
    lab.mean=lab.obs.mean.fun(dat.x.raw[dat.x.raw$patient_num%in%patient_num.deceased,], dat.calendar[dat.calendar$patient_num%in%patient_num.deceased,], nm.patient_num, nm.days_since_admission) 
    lab.mean.std=lab.obs.mean.std.fun(dat.x.raw[dat.x.raw$patient_num%in%patient_num.deceased,], dat.calendar[dat.calendar$patient_num%in%patient_num.deceased,], nm.patient_num, nm.days_since_admission) 
  }  
  res[[nm.event]]=list(lab.obs=lab.obs, lab.mean=lab.mean, lab.mean.std=lab.mean.std)
}
res
}
