data_analysis_clean=function (dir.input, code.dict, nm.event, dat.surv.raw, dat.x.raw, dat.dem.raw, 
          nm.patient_num, nm.days_since_admission, nm.value, patient.keep) 
{
  code.dict = apply(code.dict, 2, as.character)
  combine.set = c("48065-7", "48066-5")
  combine.nm = paste(combine.set, collapse = ":")
  code.dict = data.frame(rbind(code.dict, c(combine.nm, "D-dimer")))
  
  dat.event=dat.dem.raw[,c("patient_num", "admission_date", "days_since_admission","last_discharge_date", "severe_date", "death_date")]
  dat.event$x1=as.numeric(as.Date(dat.event$last_discharge_date)-as.Date(dat.event$admission_date))
  dat.event$x2=as.numeric(as.Date(dat.event$severe_date)-as.Date(dat.event$admission_date))
  dat.event$x3=as.numeric(as.Date(dat.event$death_date)-as.Date(dat.event$admission_date))
  dat.event$x1[which(dat.event$x1<0|is.na(dat.event$x1))]=dat.event$admission_date[which(dat.event$x1<0|is.na(dat.event$x1))]
  dat.event[which(dat.event<0)]=NA
  dat.event$lastevent_time=pmax(dat.event$x1, dat.event$x2, dat.event$x3, na.rm=T)
  dat.event=dat.event[,c("patient_num", "lastevent_time")]
  dat.censor=data.frame(dat.x.raw[,c("patient_num", "days_since_admission")])
  dat.censor=dat.censor[duplicated(dat.censor)!=1,]
  dat.censor=dat.censor[order(dat.censor$patient_num, dat.censor$days_since_admission),]
  dat.censor=data.table(dat.censor)
  dat.censor=data.frame(dat.censor[,lapply(.SD, max), by="patient_num"])
  colnames(dat.censor)[2]="lastcode_time"
  
  dat.censor=left_join(dat.event,dat.censor, by="patient_num")
  dat.censor$last_info=pmax(dat.censor$lastevent_time, dat.censor$lastcode_time)
  dat.censor=dat.censor[,c("patient_num", "last_info")]
  dat.surv = data_surv_clean(dat.surv.raw, nm.patient_num, 
                             nm.days_since_admission, dat.censor, daymax = 30, nm.event)
  dat.lab = data_lab_clean(dat.x.raw, code.dict, nm.patient_num, 
                           nm.days_since_admission, nm.value, day = 0)
  dat.lab=lab_outlier_fun(dat.lab)
  dat.lab=dat.lab$dat.lab.new
  dat.dem = data_dem_clean(dat.dem.raw, nm.patient_num, nm.gender = "sex", 
                           nm.age = "age_group", nm.race = "race")
  dat.analysis = left_join(dat.surv, dat.dem, by = nm.patient_num)
  dat.analysis = left_join(dat.analysis, dat.lab, by = nm.patient_num)
  dat.analysis$sex[dat.analysis$sex=="other"]="female"
  dat.analysis$age_group_new = dat.analysis$age_group
  dat.analysis$age_group_new[dat.analysis$age_group %in% c("00to02", 
                                                           "03to05", "06to11", "12to17", "18to25")] = "00to25"
  dat.analysis$age_group_new[dat.analysis$age_group %in% "other"]="50to69" ## combine other with 50to69
  
  dat.analysis$age_group_new <- factor(dat.analysis$age_group_new, 
                                       levels = c("00to25", "26to49", "50to69", "70to79", "80plus"))
  
  dat.analysis$race_new = dat.analysis$race
  dat.analysis$race_new[dat.analysis$race %in% c("American_indian", 
                                                 "Hawaiian_pacific_islander", "Hispanic_latino", "Other")] = "Hispanic and Other"
  dat.analysis$race_new <- factor(dat.analysis$race_new, levels = c("White", 
                                                                    "Black", "Asian", "Hispanic and Other"))
  dat.analysis = dat.analysis[dat.analysis$patient_num %in% 
                                patient.keep, ]
  dat.analysis
}