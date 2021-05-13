
dat.prep.fun=function(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls){
  dat=dat.survival[[paste0("dat.analysis.",nm.event)]]
  dat=dat[,duplicated(colnames(dat))!=1]
  dat.calendar=dat.survival[["dat.calendar"]]
  dat=left_join(dat.calendar[,c("patient_num", "calendar_day", "calendar_date")], dat, by="patient_num")
  nm.lab.keep.new=colnames(dat)[colnames(dat)%in%nm.lab.keep]
  dat=dat[,c("patient_num", "days_since_admission", nm.event, "calendar_date","calendar_day", nm.dem, nm.lab.keep.new, nm.cls)]
  dat$age_group_new = relevel(dat$age_group_new, ref = "50to69")
  dat$sex=tolower(dat$sex)
  dat$sex=factor(dat$sex, level=c("male", "female"))
  dat$race_new=factor(dat$race_new, level=c("White","Black","Asian","Hispanic and Other"))
  dat.lab.impute=do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];mm=mean(x,na.rm=T);x[is.na(x)]=mm; x}))
  obs.indx=1*do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];is.na(x)!=1}))
  colnames(obs.indx)=paste0("obs_", nm.lab.keep.new)
  dat[,nm.lab.keep.new]=dat.lab.impute
  dat.cls.impute=dat[,nm.cls]
  dat.cls.impute[is.na(dat.cls.impute)]=mean(dat.cls.impute, na.rm=T)
  cls.indx=1*(is.na(dat[,nm.cls])!=1)
  dat[,nm.cls]=dat.cls.impute

  dat=cbind(dat, obs.indx, obs_charlson_score=cls.indx)
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
  nm.lab.new=setdiff(colnames(dat),c("patient_num", "days_since_admission", "severedeceased", "deceased","calendar_date","calendar_day",        
                                     "age","sex","race"))
  nm.lab.new=setdiff(nm.lab.new, c(nm.lab.new[grepl("obs_",nm.lab.new)],"charlson_score"))
  #if("DD"%in%nm.lab.new){
  #dat$DD=log(dat$DD+1)}
  dat[,colnames(dat)%in%nm.lab.new]=log(dat[,colnames(dat)%in%nm.lab.new]+1)
  dat
  }

