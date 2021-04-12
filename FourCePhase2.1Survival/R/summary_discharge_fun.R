summary_discharge_fun=function(dir.input){
dat.cc=read.csv(paste0(dir.input, "LocalPatientClinicalCourse.csv"))
dat.summary=read.csv(paste0(dir.input, "LocalPatientSummary.csv"))
sum(dat.summary$last_discharge_date<dat.summary$death_date)
xx=dat.summary[which(dat.summary$last_discharge_date<dat.summary$death_date & dat.summary$last_discharge_date!="1900-01-01"),c("patient_num", "last_discharge_date", "death_date")]

dat.cc=dat.cc[order(dat.cc$days_since_admission),]
patient_num.list=sort(unique(dat.cc$patient_num))
dat.res=NULL
for(patient_num in patient_num.list){
  
  dat.cc.sub=dat.cc[dat.cc$patient_num==patient_num ,]
  dat.cc.change=dat.cc.sub %>%mutate(in_hospital.change=cumsum(c(1, diff(in_hospital) != 0)))
  
  dat.cc.change.admission=dat.cc.change[which(dat.cc.change$in_hospital==1),]
  dat.cc.change.discharge=dat.cc.change[which(dat.cc.change$in_hospital==0),]
  index.admission=which(!duplicated(dat.cc.change.admission$in_hospital.change))
  index.discharge=which(!duplicated(dat.cc.change.discharge$in_hospital.change))
  last_discharge_day=ifelse(length(index.discharge)==0,NA,dat.cc.change.discharge[index.discharge[length(index.discharge)],"days_since_admission"])-1
  last_discharge_date=as.character(as.Date(ifelse(length(index.discharge)==0, NA,dat.cc.change.discharge[index.discharge[length(index.discharge)],"calendar_date"]))-1)
  
  if((length(index.admission)-length(index.discharge))==1){index.discharge=c(index.discharge, 10000)}
  admission_discharge_day=apply(cbind(dat.cc.change.admission[index.admission,"days_since_admission"], 
                                      dat.cc.change.discharge[index.discharge,"days_since_admission"]),1, paste,collapse="->")
  
  admission_discharge_date=apply(cbind(dat.cc.change.admission[index.admission,"calendar_date"], 
                                       dat.cc.change.discharge[index.discharge,"calendar_date"]),1, paste,collapse="->")
  
  last_obs_day=dat.cc.change[dim(dat.cc.change)[1], "days_since_admission"]
  last_obs_date=dat.cc.change[dim(dat.cc.change)[1], "calendar_date"]
  
  first_obs_day=dat.cc.change[1, "days_since_admission"]
  first_obs_date=dat.cc.change[1, "calendar_date"]
  
  deseased_day=dat.cc.change[which(dat.cc.change$deceased==1)[1], "days_since_admission"]
  deseased_date=dat.cc.change[which(dat.cc.change$deceased==1)[1], "calendar_date"]
  
  dat.after=dat.cc.change[dat.cc.change$days_since_admission>last_discharge_day,]
  severe_change_after_last_discharge=ifelse(length(table(dat.after$severe))>1,1,0)
  deceased_change_after_last_discharge=ifelse(length(table(dat.after$deceased))>1,1,0)
  
  dat.tmp=data.frame(siteid=siteid, patient_num=patient_num, 
                     admission_discharge_day=admission_discharge_day, admission_discharge_date=admission_discharge_date, 
                     last_discharge_day=last_discharge_day, last_discharge_date=last_discharge_date, 
                     first_obs_day=first_obs_day,  first_obs_date=first_obs_date, 
                     last_obs_day=last_obs_day, last_obs_date=last_obs_date, 
                     deseased_day=deseased_day, deseased_date=deseased_date,
                     severe_change_after_last_discharge=severe_change_after_last_discharge, deceased_change_after_last_discharge=deceased_change_after_last_discharge) 
  dat.res=rbind(dat.res, dat.tmp)
}

dat.res=dat.res[is.na(dat.res$patient_num)!=1,]
### how many patients only have one admission
junk=table(table(dat.res$patient_num))
cat("distribution of number of admissions: \n")
print(junk)

#which(table(dat.res$patient_num)==11)
#dat.res[dat.res$patient_num==849,]
dat7=dat.res[which(is.na(dat.res$last_discharge_day)!=1 & dat.res$last_discharge_day<7),]
dat14=dat.res[which(is.na(dat.res$last_discharge_day)!=1 & dat.res$last_discharge_day<14),]
dat28=dat.res[which(is.na(dat.res$last_discharge_day)!=1 & dat.res$last_discharge_day<28),]

cat("All patients have information after discharge until the current last calendar date, even if deceased. \n")

cat("number of patients last discharge within 7 days: ", dim(dat7)[1], "\n")
cat("number of patients last discharge within 7 days and change severe status after that: ", dim(dat7[dat7$severe_change_after_last_discharge==1,])[1], "\n")
cat("number of patients last discharge within 7 days and change deceased status after that: ", dim(dat7[dat7$deceased_change_after_last_discharge==1,])[1], "\n")

cat("number of patients last discharge within 14 days: ", dim(dat14)[1], "\n")
cat("number of patients last discharge within 14 days and change severe status after that: ", dim(dat14[dat14$severe_change_after_last_discharge==1,])[1], "\n")
cat("number of patients last discharge within 14 days and change deceased status after that: ", dim(dat14[dat14$deceased_change_after_last_discharge==1,])[1], "\n")

cat("number of patients last discharge within 28 days: ", dim(dat28)[1], "\n")
cat("number of patients last discharge within 28 days and change severe status after that: ", dim(dat28[dat28$severe_change_after_last_discharge==1,])[1], "\n")
cat("number of patients last discharge within 28 days and change deceased status after that: ", dim(dat28[dat28$deceased_change_after_last_discharge==1,])[1], "\n")

cat("number of patients deceased after last discharge: ", table(dat.res$deceased_change_after_last_discharge)[2], "\n")

res=left_join(dat.res, dat.summary, by="patient_num")
#table(as.Date(xx$last_discharge_date.x)-as.Date(xx$last_discharge_date.y))
## discharge date should be the next day of date hospitalization become zero
#table(as.Date(xx$deseased_date)-as.Date(xx$death_date))
res
}
