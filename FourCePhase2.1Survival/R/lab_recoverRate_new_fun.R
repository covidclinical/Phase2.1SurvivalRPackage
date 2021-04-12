lab_recoverRate_new_fun=function(dir.input, pat.days.cut1=0, pat.days.cut2=14, days.range=c(1:14), code.dict,LocalPatientObservations, LocalPatientClinicalCourse, dat.survival,calendar.date.cut, myscale, is.ns=0){
dat.calendar=dat.survival$dat.calendar
dat.x.raw = LocalPatientObservations
dat.cc=LocalPatientClinicalCourse
dat.cc=dat.cc[dat.cc$in_hospital==1,]
dat.cc.stay=dat.cc[,c("patient_num", "days_since_admission")]
dat.cc.stay=data.frame(data.table(dat.cc.stay)[,max(days_since_admission), by="patient_num"])
colnames(dat.cc.stay)[2]="max_day"
dat.x.raw=left_join(dat.x.raw,dat.cc.stay, by="patient_num")
dat.x.raw=dat.x.raw[-which(dat.x.raw$days_since_admission>dat.x.raw$max_day),]
if(pat.days.cut2!=Inf){
  patient_num.keep=unique(dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1 & 
                                                                                 dat.survival$dat.analysis.deceased$days_since_admission>=pat.days.cut1 &
                                                                                 dat.survival$dat.analysis.deceased$days_since_admission<=pat.days.cut2)])
}else{
  patient_num.rm=unique(dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1 & 
                                                                               dat.survival$dat.analysis.deceased$days_since_admission<pat.days.cut1)])
  patient_num.keep=unique(setdiff(dat.survival$dat.analysis.deceased$patient_num, patient_num.rm))
  
}

junk=NULL
for(myday in days.range){
dat.lab.tmp=data_lab_clean2(code.dict, dat.x.raw, 
                            nm.value="value", day=myday) 
junk[[myday]]=dat.lab.tmp
}

#nm.lab.all=colnames(junk[[1]])[-1]
nm.lab.all=c("total_bilirubin","creatinine", "Ferritin", "D_dimer", "C_reactive_protein_CRP_Normal_Sensitivity", "albumin","neutrophil_count")
nm.lab.all=setdiff(nm.lab.all,"cardiac_troponin_Normal_Sensitivity")
nm.lab.all=nm.lab.all[nm.lab.all%in%unique(unlist(lapply(days.range, function(myday) colnames(junk[[myday]]))))]

dat.fit=NULL
for(nm.lab in nm.lab.all){
  tmp=NULL
for(myday in days.range){
  if(nm.lab%in%colnames(junk[[myday]])){
  if(myscale=="original"){tmp=rbind(tmp,data.frame(day=myday,patient_num=junk[[myday]][,"patient_num"], y=junk[[myday]][,nm.lab]))}
  if(myscale=="log"){tmp=rbind(tmp,data.frame(day=myday,patient_num=junk[[myday]][,"patient_num"], y=log(junk[[myday]][,nm.lab])))}
  }else{
    if(myscale=="original"){tmp=rbind(tmp,data.frame(day=myday,patient_num=junk[[myday]][,"patient_num"], y=NA))}
    if(myscale=="log"){tmp=rbind(tmp,data.frame(day=myday,patient_num=junk[[myday]][,"patient_num"], y=NA))}
  }
  }
tmp=left_join(tmp, dat.calendar,by="patient_num")
dat.fit[[nm.lab]]=tmp
}

res=NULL
for(nm.lab in nm.lab.all){
  print(nm.lab)
  dat.lab=dat.fit[[nm.lab]]
  dat.lab=dat.lab[complete.cases(dat.lab),]
  dat.lab=dat.lab[dat.lab$y%in%c(Inf, -Inf)!=1,]
  dat.lab=dat.lab[dat.lab$patient_num%in%patient_num.keep,]
 if(dim(dat.lab)[1]>10){
  if(is.ns==0){
  m1.early <- tryCatch(lmer(y ~ day+(day|patient_num), data=dat.lab[dat.lab$calendar_date<calendar.date.cut,]),error=function(e) NA)
  m1.late <- tryCatch(lmer(y ~ day+(day|patient_num), data=dat.lab[dat.lab$calendar_date>=calendar.date.cut,]),error=function(e) NA)

  s1.early=NA
  s1.late=NA
  }else{
  myknots=c(min(dat.lab$day)+2, floor(mean(range(dat.lab$day))), max(dat.lab$day)-3)
  m1.early<- tryCatch(lmer(y ~ ns(day,knots=myknots)+(day|patient_num), data=dat.lab[dat.lab$calendar_date<calendar.date.cut,]),error=function(e) NA)
  m1.late <- tryCatch(lmer(y ~ ns(day,knots=myknots)+(day|patient_num), data=dat.lab[dat.lab$calendar_date>=calendar.date.cut,]),error=function(e) NA)

  ff=ns(c(days.range,dat.lab[dat.lab$calendar_date<calendar.date.cut,]$day),knots=myknots)[1:length(days.range),]

  s1.early=tryCatch(cbind(1,ff)%*%coef(summary(m1.early))[, "Estimate"],error=function(e) NA)
  s1.late=tryCatch(cbind(1,ff)%*%coef(summary(m1.late))[, "Estimate"], error=function(e) NA)
  }
 # plot: a1f1+a2f2+a3f3
 # a1,a2,a3 covariance matrix
 # knots: 7,14,21; cubic spline
 # day up to 30 days
 # risk prediction
 # 
  ss.early=tryCatch(sd(ranef(m1.early)$patient_num[,2]),error=function(e) NA)
  mm.early=tryCatch(coef(summary(m1.early))[-1, "Estimate"],error=function(e) NA)
  vcov.early=tryCatch((summary(m1.early))$vcov, error=function(e) NA)

  ss.late=tryCatch(sd(ranef(m1.late)$patient_num[,2]),error=function(e) NA)
  mm.late=tryCatch(coef(summary(m1.late))[-1, "Estimate"],error=function(e) NA)
  vcov.late=tryCatch((summary(m1.late))$vcov, error=function(e) NA)
  
  randomeffect.early=tryCatch(hist(ranef(m1.early)$patient_num[,2], breaks=100, freq=0), error=function(e) NA)
  randomeffect.late=tryCatch(hist(ranef(m1.late)$patient_num[,2], breaks=100, freq=0),error=function(e) NA)
  
  res[[nm.lab]]=list(
    ss.early=ss.early, 
    mm.early=mm.early, 
    vcov.early=vcov.early,
    s1.early=s1.early,
    ss.late=ss.late, 
    mm.late=mm.late, 
    vcov.late=vcov.late,
    s1.late=s1.late#,
    #randomeffect.early=randomeffect.early, 
    #randomeffect.late=randomeffect.late
    )
 }
}
res
}
