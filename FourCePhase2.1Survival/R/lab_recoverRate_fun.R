lab_recoverRate_fun=function(dir.input, code.dict,LocalPatientObservations, dat.survival,calendar.date.cut, myscale, is.ns=0){
dat.calendar=dat.survival$dat.calendar
dat.x.raw = LocalPatientObservations
patient_num.severe=dat.survival$dat.analysis.severe$patient_num[which(dat.survival$dat.analysis.severe$severe==1)]
patient_num.deceased=dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1)]
patient_num.severedeceased=dat.survival$dat.analysis.severedeceased$patient_num[which(dat.survival$dat.analysis.severedeceased$severedeceased==1)]

junk=NULL
for(myday in 1:30){
dat.lab.tmp=data_lab_clean2(code.dict, dat.x.raw, 
                            nm.value="value", day=myday) 
junk[[myday]]=dat.lab.tmp
}

#nm.lab.all=colnames(junk[[1]])[-1]
nm.lab.all=c("total_bilirubin","creatinine", "Ferritin", "D_dimer", "C_reactive_protein_CRP_Normal_Sensitivity", "albumin","neutrophil_count")
nm.lab.all=setdiff(nm.lab.all,"cardiac_troponin_Normal_Sensitivity")
nm.lab.all=nm.lab.all[nm.lab.all%in%colnames(junk[[1]])]
dat.fit=NULL
for(nm.lab in nm.lab.all){
  tmp=NULL
for(myday in 1:30){
  if(myscale=="original"){tmp=rbind(tmp,data.frame(day=myday,patient_num=junk[[myday]][,"patient_num"], y=junk[[myday]][,nm.lab]))}
  if(myscale=="log"){tmp=rbind(tmp,data.frame(day=myday,patient_num=junk[[myday]][,"patient_num"], y=log(junk[[myday]][,nm.lab])))}
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
  
  if(is.ns==0){
  m1.early <- tryCatch(lmer(y ~ day+(day|patient_num), data=dat.lab[dat.lab$calendar_date<calendar.date.cut,]),error=function(e) NA)
  m1.late <- tryCatch(lmer(y ~ day+(day|patient_num), data=dat.lab[dat.lab$calendar_date>=calendar.date.cut,]),error=function(e) NA)
  s1.early=NA
  s1.late=NA
  }else{
  m1.early <- tryCatch(lmer(y ~ ns(day,knots=c(7,14,21))+(day|patient_num), data=dat.lab[dat.lab$calendar_date<calendar.date.cut,]),error=function(e) NA)
  m1.late <- tryCatch(lmer(y ~ ns(day,knots=c(7,14,21))+(day|patient_num), data=dat.lab[dat.lab$calendar_date>=calendar.date.cut,]),error=function(e) NA)
  ff=ns(c(0:30,dat.lab[dat.lab$calendar_date<calendar.date.cut,]$day),knots=c(7,14,21))[0:30,]

  s1.early=tryCatch(ff%*%coef(summary(m1.early))[-1, "Estimate"],error=function(e) NA)
  s1.late=tryCatch(ff%*%coef(summary(m1.late))[-1, "Estimate"], error=function(e) NA)
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
  
  #randomeffect.early=tryCatch(hist(ranef(m1.early)$patient_num[,2], breaks=100), error=function(e) NA)
  #randomeffect.late=tryCatch(hist(ranef(m1.late)$patient_num[,2], breaks=100),error=function(e) NA)
  
  res[[nm.lab]]=list(ss.early=ss.early, mm.early=mm.early, vcov.early=vcov.early,s1.early=s1.early,
                     ss.late=ss.late, mm.late=mm.late, vcov.late=vcov.late,s1.late=s1.late#,
                     #randomeffect.early=randomeffect.early, 
                     #randomeffect.late=randomeffect.late
                     )
}
res
}
