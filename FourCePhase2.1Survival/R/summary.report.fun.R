summary.report.fun=function(dat.survival, siteid, dir.output){
  dat.calendar=dat.survival[["dat.calendar"]]
  dat.severe=dat.survival[["dat.analysis.severe"]]
  dat.severedeceased=dat.survival[["dat.analysis.severedeceased"]]
  dat.deceased=dat.survival[["dat.analysis.deceased"]]
  
  tab.all=table(substr(dat.calendar$calendar_date,1,7))
    
  cat("Calendar Date:", as.character(sort(dat.calendar$calendar_date)[1]), " to ", as.character(sort(dat.calendar$calendar_date)[dim(dat.calendar)[1]]), "\n")
  cat("Number of patients:", dim(dat.calendar)[1], "\n")
  cat("Number of severe patients:", sum(dat.severe$severe), "\n")
  cat("Number of severe patients at day0:", sum(dat.severe$severe==1 & dat.severe$days_since_admission==0), "\n")
  cat("Number of severe patients admitted on different calendar months:\n")
  tmp=left_join(dat.calendar,dat.severe,by="patient_num")
  tab.severe=table(substr(tmp$calendar_date,1,7), tmp$severe)
  colnames(tab.severe)=paste0("severe",colnames(tab.severe))
  print(tab.severe)
  
  cat("Number of dead patients:", sum(dat.deceased$deceased), "\n")
  cat("Number of dead patients at day0:", sum(dat.deceased$deceased==1 & dat.deceased$days_since_admission==0), "\n")
  cat("Number of dead patients admitted on different calendar months:\n")
  tmp=left_join(dat.calendar,dat.deceased,by="patient_num")
  tab.deceased=table(substr(tmp$calendar_date,1,7), tmp$deceased)
  colnames(tab.deceased)=paste0("deceased",colnames(tab.deceased))
  print(tab.deceased)

  tmp=left_join(dat.calendar,dat.severe[,c("patient_num", "severe")],by="patient_num")
  tmp=left_join(tmp, dat.deceased[,c("patient_num","deceased")], by="patient_num")
  tab.severedeceased=table(substr(tmp$calendar_date,1,7), 1*((tmp$severe+tmp$deceased)>0))
  colnames(tab.severedeceased)=paste0("severedeceased",colnames(tab.severedeceased))
  print(tab.severedeceased)
  
  cat("Age:\n")
  print(table(dat.severe$age_group))
  cat("Sex:\n")
  print(table(dat.severe$sex))
  cat("Race:\n")
  print(table(dat.severe$race))
  res=cbind(N=as.numeric(tab.all), tab.severe, tab.deceased, tab.severedeceased)
res
}
