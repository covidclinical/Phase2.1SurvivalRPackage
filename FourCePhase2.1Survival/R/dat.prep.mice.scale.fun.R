
dat.prep.mice.scale.fun=function(dat.survival, nm.event, nm.dem, nm.lab.keep, nm.cls, method.impute, myscale,is.ind, mice.time){
  dat=dat.survival[[paste0("dat.analysis.",nm.event)]]
  dat=dat[,duplicated(colnames(dat))!=1]
  junk.rm=apply(dat,2,function(x)mean(is.na(x)))
  nm.rm=names(junk.rm)[junk.rm>0.95]
  dat=dat[,setdiff(colnames(dat), nm.rm)]
  dat.calendar=dat.survival[["dat.calendar"]]
  dat=left_join(dat.calendar[,c("patient_num", "calendar_day", "calendar_date")], dat, by="patient_num")
  nm.lab.keep.new=colnames(dat)[colnames(dat)%in%nm.lab.keep]
  dat=dat[,c("patient_num", "days_since_admission", nm.event, "calendar_date","calendar_day", nm.dem, nm.lab.keep.new, nm.cls)]
  dat$age_group_new = relevel(dat$age_group_new, ref = "50to69")
  dat$sex=tolower(dat$sex)
  dat$sex=factor(dat$sex, level=c("male", "female"))
  dat$race_new=factor(dat$race_new, level=c("White","Black","Asian","Hispanic and Other"))
  
  if(method.impute=="zero"){
    dat.lab.impute=do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];x[is.na(x)]=0; x}))
  }
  
  if(method.impute=="mean"){
  dat.lab.impute=do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];mm=mean(x,na.rm=T);x[is.na(x)]=mm; x}))
  }
  if(method.impute=="median"){
  dat.lab.impute=do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];mm=median(x,na.rm=T);x[is.na(x)]=mm; x}))
  }
  if(method.impute=="mice"){
    dat.lab.impute=do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];mm=median(x,na.rm=T);x[is.na(x)]=mm; x}))
    mice_imputes0 = mice(dat[,-c(1:5)], m=mice.time, maxit = 40, seed=1234, print=FALSE)
    mice_imputes=lapply(1:mice.time, function(ll) complete(mice_imputes0, ll))
    mice_imputes=Reduce("+", mice_imputes)/length(mice_imputes)
    dat.lab.impute= mice_imputes[,nm.lab.keep.new]
  }  
  
  mis.indx=1*do.call(cbind,lapply(nm.lab.keep.new, function(ll){x=dat[,ll];is.na(x)==1}))
  colnames(mis.indx)=paste0("mis_", nm.lab.keep.new)
  dat[,nm.lab.keep.new]=dat.lab.impute
  dat.cls.impute=dat[,nm.cls]
  if(method.impute=="zero"){
  dat.cls.impute[is.na(dat.cls.impute)]=0}
  if(method.impute=="mean"){
    dat.cls.impute[is.na(dat.cls.impute)]=mean(dat.cls.impute, na.rm=T)}
  if(method.impute=="median"){
    dat.cls.impute[is.na(dat.cls.impute)]=median(dat.cls.impute, na.rm=T)}
  if(method.impute=="mice"){
    dat.cls.impute= mice_imputes[,nm.cls]
    
    }
  cls.indx=1*(is.na(dat[,nm.cls])==1)
  dat[,nm.cls]=dat.cls.impute

  dat=cbind(dat, mis.indx, mis_charlson_score=cls.indx)
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
  nm.lab.new=setdiff(nm.lab.new, c(nm.lab.new[grepl("mis_",nm.lab.new)],"charlson_score"))
  if(myscale=="log"){
  dat[,colnames(dat)%in%nm.lab.new]=log(dat[,colnames(dat)%in%nm.lab.new]+0.5)}else{
    nm.trans=c("ALT", "AST", "CRP", "DD")
    dat[,colnames(dat)%in%nm.trans]=log(dat[,colnames(dat)%in%nm.trans]+0.5)
  }
  if(is.ind!=1){dat=dat[,colnames(dat)[grepl("mis_", colnames(dat))!=1]]}
  dat
  }

