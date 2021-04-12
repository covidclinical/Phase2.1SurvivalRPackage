data_analysis_clean=function (dir.input, code.dict, nm.event, dat.surv.raw, dat.x.raw, dat.dem.raw, 
          nm.patient_num, nm.days_since_admission, nm.value, patient.keep) 
{
  code.dict = apply(code.dict, 2, as.character)
  combine.set = c("48065-7", "48066-5")
  combine.nm = paste(combine.set, collapse = ":")
  code.dict = data.frame(rbind(code.dict, c(combine.nm, "D-dimer")))
  dat.surv = data_surv_clean(dat.surv.raw, nm.patient_num, 
                             nm.days_since_admission, daymax = 30, nm.event)
  dat.lab = data_lab_clean(dat.x.raw, code.dict, nm.patient_num, 
                           nm.days_since_admission, nm.value, day = 0)
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