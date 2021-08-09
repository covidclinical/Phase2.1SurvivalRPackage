getDCData=function (dir.input, code.dict, siteid, 
                    LocalPatientClinicalCourse, 
                    LocalPatientObservations,
                    LocalPatientSummary) 
{
  nm.patient_num = "patient_num"
  nm.days_since_admission = "days_since_admission"
  nm.event="deceased"
  nm.value = "value"
  nm.gender = "sex"
  nm.age = "age_group"
  nm.race = "race"
  dat.surv.raw =LocalPatientClinicalCourse
  dat.x.raw = LocalPatientObservations
  dat.dem.raw = LocalPatientSummary
  
  
  # define calendar_time
  dat.calendar = dat.surv.raw[dat.surv.raw$days_since_admission ==  0, c("patient_num", "calendar_date")]
  patient.keep = unique(dat.surv.raw[, "patient_num"])
  dat.calendar$calendar_time = as.Date(dat.calendar$calendar_date)
  dat.calendar$calendar_time = as.numeric(dat.calendar$calendar_time - min(dat.calendar$calendar_time) ) # as.Date("2020-03-1")
  
  
  # dat.analysis.deceased = data_analysis_clean(dir.input, code.dict, nm.event = "deceased", 
  #                                             dat.surv.raw, dat.x.raw, dat.dem.raw, nm.patient_num, 
  #                                             nm.days_since_admission, nm.value, patient.keep)
  
  # collect survival data
  # dat.surv = data_surv_clean(dat.surv.raw, nm.patient_num, nm.days_since_admission, daymax = 30, nm.event)
  daymax = 400
  dat.sub = dat.surv.raw[dat.surv.raw[, nm.days_since_admission] <= daymax, ]
  patient_num.list = sort(unique(dat.sub[, nm.patient_num]))
  data.surv.clean.kern = function(dat.sub, patient_num) {
    dat.tmp = dat.sub[dat.sub[, nm.patient_num] == patient_num, ]
    dat.event.tmp = dat.tmp[which(dat.tmp[, nm.event] ==  1), ]
    if (dim(dat.event.tmp)[1] != 0) {
      dat.end.tmp = dat.event.tmp[which.min(dat.event.tmp[,nm.days_since_admission]), c(nm.patient_num,nm.days_since_admission, nm.event)]
    }
    else {
      dat.end.tmp = dat.tmp[which.max(dat.tmp[, nm.days_since_admission]),c(nm.patient_num, nm.days_since_admission, nm.event)]
    }
    dat.end.tmp
  }
  dat.surv = do.call(rbind, lapply(patient_num.list, function(xx) data.surv.clean.kern(dat.sub,xx)))
  rownames(dat.surv) = NULL
  # save as tempdata2
  
  
  # collect lab data
  # dat.lab = data_lab_clean(dat.x.raw, code.dict, nm.patient_num, nm.days_since_admission, nm.value, day = 0)
  day = 0
  code.dict = apply(code.dict, 2, as.character)
  code.dict = data.frame(rbind(code.dict, c("48065-7:48066-5", "D-dimer")))
  
  dat.sub = dat.x.raw[dat.x.raw[, nm.days_since_admission] == day & dat.x.raw[, "concept_type"] %in% c("LAB-LOINC"), ]
  dat.sub$concept = paste(dat.sub$concept_type, dat.sub$concept_code, sep = ":")
  dat.sub = dat.sub[, c(nm.patient_num, "concept", nm.value)]
  
  dat.sub.wide <- dat.sub %>%pivot_wider(values_from = value, names_from = concept, values_fn = mean)
  dat.sub.wide[dat.sub.wide == -99] = NA
  dat.sub.wide[dat.sub.wide == -999] = NA
  if ("LAB-LOINC:48065-7" %in% colnames(dat.sub.wide)) {
    dat.sub.wide[, "LAB-LOINC:48065-7:48066-5"] = 0.5 * dat.sub.wide[, "LAB-LOINC:48065-7"]
  }
  if ("LAB-LOINC:48066-5" %in% colnames(dat.sub.wide)) {
    dat.sub.wide[, "LAB-LOINC:48065-7:48066-5"] = dat.sub.wide[,  "LAB-LOINC:48066-5"]
  }
  dat.sub.wide = dat.sub.wide[, setdiff(colnames(dat.sub.wide),c("LAB-LOINC:48065-7", "LAB-LOINC:48066-5"))] # 1+20
  
  nm.lab = colnames(dat.sub.wide)[-1]
  code.dict.new = code.dict
  colnames(code.dict.new)[1] = "concept"
  code.dict.new[, 1] = paste("LAB-LOINC", code.dict.new[, 1], sep = ":")
  # code.dict.new = apply(code.dict.new, 2, as.character)
  # code.dict.new = data.frame(rbind(code.dict.new, cbind(nm.lab[grepl("MED-", nm.lab)==1],nm.lab[grepl("MED-", nm.lab)==1]) ) )
  nm.lab = left_join(data.frame(concept = nm.lab), code.dict.new,by = "concept")
  nm.lab = nm.lab[, 2]
  nm.lab = gsub(" ", "_", nm.lab); nm.lab = gsub("\\(", "", nm.lab); nm.lab = gsub("\\)", "", nm.lab); nm.lab = gsub("-", "_", nm.lab)
  colnames(dat.sub.wide) = c(colnames(dat.sub.wide)[1], nm.lab)
  rownames(dat.sub.wide) = NULL
  dat.lab=data.frame(dat.sub.wide[, which(is.na(colnames(dat.sub.wide)) != 1)])
  
  # collect med data
  dat.sub = dat.x.raw[dat.x.raw[, "concept_type"] %in% c("MED-CLASS"), ]
  dat.sub$concept = paste(dat.sub$concept_type, dat.sub$concept_code, sep = ":")
  dat.sub = dat.sub[, c(nm.patient_num, "concept","value")]
  dat.sub.wide <- dat.sub %>%pivot_wider(values_from = value, names_from = concept)
  for (i in 2:ncol(dat.sub.wide)){
    dat.sub.wide[[i]]=sapply(1:nrow(dat.sub.wide), function(kk){ length((dat.sub.wide[[i]])[[kk]]) } )
  }
  dat.med=data.frame(dat.sub.wide)
  
  
  # collect demographic data
  # dat.dem = data_dem_clean(dat.dem.raw, nm.patient_num, nm.gender = "sex", nm.age = "age_group", nm.race = "race")
  dat.sub = data.frame(dat.dem.raw[, c(nm.patient_num, nm.gender, nm.age, nm.race)])
  dat.dem=dat.sub
  
  
  # join three datasets
  dat.analysis = left_join(dat.surv, dat.dem, by = nm.patient_num)
  dat.analysis = left_join(dat.analysis, dat.lab, by = nm.patient_num)
  dat.analysis = left_join(dat.analysis, dat.med, by = nm.patient_num)
  
  dat.analysis$sex[dat.analysis$sex=="other"]="female"
  dat.analysis$age_group_new = dat.analysis$age_group
  dat.analysis$age_group_new[dat.analysis$age_group %in% c("00to02", "03to05", "06to11", "12to17", "18to25")] = "00to25"
  dat.analysis$age_group_new[dat.analysis$age_group %in% "other"]="50to69" ## combine other with 50to69
  dat.analysis$age_group_new <- factor(dat.analysis$age_group_new, levels = c("00to25", "26to49", "50to69", "70to79", "80plus"))
  
  dat.analysis$race_new = dat.analysis$race
  dat.analysis$race_new[dat.analysis$race %in% c("american_indian", "hawaiian_pacific_islander", "hispanic_latino", "other")] = "hispanic and other"
  dat.analysis$race_new <- factor(dat.analysis$race_new, levels = c("white",  "black", "asian", "hispanic and other"))
  
  dat.analysis = dat.analysis[dat.analysis$patient_num %in%  patient.keep, ]
  dat.analysis.deceased=dat.analysis
  dat.survival0 = list(dat.calendar = dat.calendar, dat.analysis.deceased = dat.analysis.deceased)
  # save as tempdata3
  
  
  # add comorbidity score
  comorb=map_charlson_codes(LocalPatientObservations)
  index_scores <- comorb[[3]]
  dat.cls0<- data.frame(index_scores %>% dplyr::select(patient_num, charlson_score))
  dat.cls0$patient_num=as.numeric(dat.cls0$patient_num)
  
  dat.survival=dat.survival0
  dat.cls=dat.cls0
  dat.survival$dat.analysis.deceased=left_join(dat.survival$dat.analysis.deceased, dat.cls, by="patient_num")
  
  
  dat.survival$dat.calendar$calendar_date=substr(dat.survival$dat.calendar$calendar_date,1,7)
  dat.survival$dat.calendar=dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date<="2021-01",] #"2021-03"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date<"2020-03"]="2020-03"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-04"]="2020-03"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-06"]="2020-05"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-08"]="2020-07"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date=="2020-10"]="2020-09"
  dat.survival$dat.calendar$calendar_date[dat.survival$dat.calendar$calendar_date>="2020-12"]="2020-11"
  
  # patient_num.out=unique(dat.survival$dat.calendar[dat.survival$dat.calendar$calendar_date>"2020-10","patient_num"])
  # dat.survival=rmOutlierSurvivalData(dat.survival,patient_num.out)
  
  
  #### further clean data
  # setting up variable names
  nm.dem.new=c("sex", "age_group_new", "race_new")
  nm.dem=c("sex", "age", "race")
  nm.cls="charlson_score"
  
  nm.lab.LabAll=colnames(dat.survival$dat.analysis.deceased)
  nm.lab.LabAll=nm.lab.LabAll[nm.lab.LabAll%in%c("age_group", "age_cat", "race_bin", nm.cls, nm.dem, nm.dem.new,
                                                 "days_since_admission", "patient_num", "value", nm.event)!=1]
  # nm.lab.LabCommon=nm.lab.LabAll
  nm.lab.LabCommon=c("alanine_aminotransferase_ALT",
                     "albumin",
                     "aspartate_aminotransferase_AST",
                     "creatinine",
                     "C_reactive_protein_CRP_Normal_Sensitivity",
                     "total_bilirubin",
                     "white_blood_cell_count_Leukocytes",
                     "lymphocyte_count",
                     "neutrophil_count","D_dimer")
  nm.lab.keep=c(nm.lab.LabCommon ) #, 
  ind=which(apply(dat.med>0,2, mean)>0.1)
  nm.med.keep=(colnames(dat.med)[ind])[-1]
  
  # period.train="all"
  # period.valid="all"
  # calendar.date.cut="2020-07"
  # t0.all=c(1:14)
  # include.ind=F   
  # include.dem=T
  # include.cls=T
  # include.lab=T
  
  # source('P:/ORD_Cho_201803054D/Wang, Xuan/dataclean/rubish/dat.prep.fun.R', echo=TRUE)
  # if(include.ind==T){dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls)}
  # if(include.ind==F){dat0=dat.prep.fun(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls); dat0=dat0[,grepl("obs_", colnames(dat0))!=1]}
  
  dat=dat.survival[[paste0("dat.analysis.",nm.event)]]
  dat=dat[,duplicated(colnames(dat))!=1]
  dat.calendar=dat.survival[["dat.calendar"]]
  dat=left_join(dat.calendar[,c("patient_num", "calendar_time", "calendar_date")], dat, by="patient_num")
  
  dat=dat[,c("patient_num", "days_since_admission", nm.event, "calendar_date","calendar_time", nm.dem.new, nm.lab.keep, nm.cls,nm.med.keep)]
  dat$age_group_new = relevel(dat$age_group_new, ref = "50to69")
  dat$sex=factor(dat$sex, level=c("male", "female"))
  dat$race_new=factor(dat$race_new, level=c("white","black","asian","hispanic and other"))
  
  dat.lab.impute=do.call(cbind,lapply(c(nm.lab.keep,nm.med.keep), function(ll){x=dat[,ll];mm=mean(x,na.rm=T);x[is.na(x)]=mm; x}))
  # obs.indx=1*do.call(cbind,lapply(nm.lab.keep, function(ll){x=dat[,ll];is.na(x)!=1}))
  # colnames(obs.indx)=paste0("obs_", nm.lab.keep)
  dat[,c(nm.lab.keep,nm.med.keep)]=dat.lab.impute
  
  dat.cls.impute=dat[,nm.cls]
  dat.cls.impute[is.na(dat.cls.impute)]=mean(dat.cls.impute, na.rm=T)
  # cls.indx=1*(is.na(dat[,nm.cls])!=1)
  dat[,nm.cls]=dat.cls.impute
  
  # dat=cbind(dat, obs.indx, obs_charlson_score=cls.indx)
  
  colnames(dat)=gsub("C_reactive_protein_CRP_Normal_Sensitivity", "CRP", colnames(dat))
  colnames(dat)=gsub("lymphocyte_count", "LYM", colnames(dat))
  colnames(dat)=gsub("D_dimer", "DD", colnames(dat))
  colnames(dat)=gsub("alanine_aminotransferase_ALT", "ALT", colnames(dat))
  colnames(dat)=gsub("aspartate_aminotransferase_AST", "AST", colnames(dat))
  colnames(dat)=gsub("total_bilirubin", "TB", colnames(dat))
  colnames(dat)=gsub("lactate_dehydrogenase_LDH", "LDH", colnames(dat))
  colnames(dat)=gsub("prothrombin_time_PT", "PT", colnames(dat))
  colnames(dat)=gsub("white_blood_cell_count_Leukocytes", "WBC", colnames(dat))
  colnames(dat)=gsub("age_group_new", "age", colnames(dat))
  colnames(dat)=gsub("race_new", "race", colnames(dat))
  colnames(dat)=gsub("cardiac_troponin_High_Sensitivity", "cardiac_high", colnames(dat))
  colnames(dat)=gsub("cardiac_troponin_Normal_Sensitivity", "cardiac_morm", colnames(dat))
  colnames(dat)=gsub("_CLASS", "", colnames(dat))
  colnames(dat)=gsub(".CLASS", "", colnames(dat))
  
  nm.lab.new=setdiff(colnames(dat),c("patient_num", "days_since_admission", "severedeceased", "deceased","calendar_date","calendar_time",        
                                     "age","sex","race"))
  nm.lab.new=setdiff(nm.lab.new, c(nm.lab.new[grepl("obs_",nm.lab.new)],"charlson_score"))
  dat[,colnames(dat)%in%nm.lab.new]=log(dat[,colnames(dat)%in%nm.lab.new]+1)
  dat0=dat
  # dat0=dat0[,grepl("obs_", colnames(dat0))!=1]
  
  dat0=dat0[dat0$days_since_admission!=0,]
  dat0
}
