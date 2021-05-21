dem.report.fun=function(dat.survival, siteid, dir.output){
  dat.calendar=dat.survival[["dat.calendar"]]
  dat.dem=dat.survival[["dat.analysis.severe"]]
  dat.dem=dat.dem[,c("patient_num", "sex", "age_group_new", "race_new")]
  tmp=left_join(dat.calendar,dat.dem,by="patient_num")
  tmp.sex=prop.table(table(tmp[,c("calendar_date", "sex")]),1)
  tmp.age=prop.table(table(tmp[,c("calendar_date", "age_group_new")]),1)
  tmp.race=prop.table(table(tmp[,c("calendar_date", "race_new")]),1)
res=list(res.sex=tmp.sex, res.age=tmp.age, res.race=tmp.race)
res
}
