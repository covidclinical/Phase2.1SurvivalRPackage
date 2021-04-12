
dat.prep.outcome.R1.fun=function(dat.survival, nm.event){
  dat=dat.survival[[paste0("dat.analysis.",nm.event)]]
  dat=dat[,duplicated(colnames(dat))!=1]
  dat.calendar=dat.survival[["dat.calendar"]]
  dat=left_join(dat.calendar[,c("patient_num", "calendar_day", "calendar_date")], dat, by="patient_num")
  dat=dat[,c("patient_num", "days_since_admission", nm.event, "calendar_date","calendar_day")]
  dat
  }

