getSurvivalData=function (dir.input, code.dict, siteid) 
{
  nm.patient_num = "patient_num"
  nm.days_since_admission = "days_since_admission"
  nm.value = "value"
  dat.surv.raw = read.csv(paste0(dir.input, "/LocalPatientClinicalCourse.csv"))
  dat.surv.raw$severedeceased = ifelse((dat.surv.raw$severe + 
                                          dat.surv.raw$deceased) > 0, 1, 0)
  dat.x.raw = read.csv(paste0(dir.input, "/LocalPatientObservations.csv"))
  dat.dem.raw = read.csv(paste0(dir.input, "/LocalPatientSummary.csv"))
  dat.calendar = dat.surv.raw[dat.surv.raw$days_since_admission == 
                                0, c("patient_num", "calendar_date")]
  patient.keep = unique(dat.surv.raw[, "patient_num"])
  dat.calendar$calendar_day = as.Date(dat.calendar$calendar_date)
  dat.calendar$calendar_day = as.numeric(dat.calendar$calendar_day - 
                                           min(dat.calendar$calendar_day))
  dat.calendar$calendar_date <- as.character(dat.calendar$calendar_date)
  
  
  dat.analysis.severe = data_analysis_clean(dir.input, code.dict, nm.event = "severe", 
                                            dat.surv.raw, dat.x.raw, dat.dem.raw, nm.patient_num, 
                                            nm.days_since_admission, nm.value, patient.keep)
  dat.analysis.deceased = data_analysis_clean(dir.input, code.dict, nm.event = "deceased", 
                                              dat.surv.raw, dat.x.raw, dat.dem.raw, nm.patient_num, 
                                              nm.days_since_admission, nm.value, patient.keep)
  dat.analysis.severedeceased = data_analysis_clean(dir.input, code.dict,
                                                    nm.event = "severedeceased", dat.surv.raw, dat.x.raw, 
                                                    dat.dem.raw, nm.patient_num, nm.days_since_admission, 
                                                    nm.value, patient.keep)
  dat.lab = read.csv(paste0(dir.input, "/Labs-", siteid, ".csv"))
  dat.med = read.csv(paste0(dir.input, "/Medications-", siteid, 
                            ".csv"))
  dat.diag = read.csv(paste0(dir.input, "/Diagnoses-", siteid, 
                             ".csv"))
  dat.dem = read.csv(paste0(dir.input, "/Demographics-", siteid, 
                            ".csv"))
  dat.dc = read.csv(paste0(dir.input, "/DailyCounts-", siteid, 
                           ".csv"))
  dat.cc = read.csv(paste0(dir.input, "/ClinicalCourse-", siteid, 
                           ".csv"))
  dat.survival = list(dat.calendar = dat.calendar, dat.analysis.severe = dat.analysis.severe, 
                      dat.analysis.deceased = dat.analysis.deceased, dat.analysis.severedeceased = dat.analysis.severedeceased, 
                      dat.lab = dat.lab, dat.med = dat.med, dat.diag = dat.diag, 
                      dat.dem = dat.dem, dat.dc = dat.dc, dat.cc = dat.cc)
  dat.survival
}
