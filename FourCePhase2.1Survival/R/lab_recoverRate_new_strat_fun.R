lab_recoverRate_new_strat_fun=function(dir.input, strat="max_day", dat.cls, pat.days.cut1=0, pat.days.cut2=14, days.range=c(1:14), code.dict, nm.lab.keep,LocalPatientObservations,LocalPatientClinicalCourse, dat.survival,calendar.date.cut, myscale, is.ns=0){
dat.calendar=dat.survival$dat.calendar
dat.x.raw = LocalPatientObservations
dat.cc=LocalPatientClinicalCourse
dat.cc=dat.cc[dat.cc$in_hospital==1,]
dat.cc.stay=dat.cc[,c("patient_num", "days_since_admission")]
dat.cc.stay=data.frame(data.table(dat.cc.stay)[,max(days_since_admission), by="patient_num"])
colnames(dat.cc.stay)[2]="max_day"
summary.stay=hist(dat.cc.stay$max_day,plot=F)

dat.x.raw=left_join(dat.x.raw,dat.cc.stay, by="patient_num")
dat.x.raw=dat.x.raw[-which(dat.x.raw$days_since_admission>dat.x.raw$max_day),]
if(pat.days.cut2!=Inf){
patient_num.keep=unique(dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1 & 
                                                                               dat.survival$dat.analysis.deceased$days_since_admission>=pat.days.cut1 &
                                                                               dat.survival$dat.analysis.deceased$days_since_admission<=pat.days.cut2)])
}else{
patient_num.rm=unique(dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1 & 
                                                                                 dat.survival$dat.analysis.deceased$days_since_admission<=pat.days.cut1)])
patient_num.keep=unique(setdiff(dat.survival$dat.analysis.deceased$patient_num, patient_num.rm))
  
}

#### patient_num summary
dat.calendar.keep=dat.calendar[dat.calendar$patient_num%in%patient_num.keep,]
day.range.list=list()
day.range.list[["week1"]]=c(0:7)
day.range.list[["week2"]]=c(8:14)
day.range.list[["week3"]]=c(15:60)
day.range.list[["all"]]=c(0:60)

dat.calendar.early=dat.calendar[dat.calendar$calendar_date<calendar.date.cut,]
dat.calendar.late=dat.calendar[dat.calendar$calendar_date>=calendar.date.cut,]
dat.calendar.early=left_join(dat.calendar.early,dat.cc.stay, by="patient_num")
dat.calendar.late=left_join(dat.calendar.late,dat.cc.stay, by="patient_num")

resN=NULL
for(ii in names(day.range.list)){
  myrange=day.range.list[[ii]]
  n.early=dim(dat.calendar.early[dat.calendar.early$max_day%in%myrange,])[1]
  n.late=dim(dat.calendar.late[dat.calendar.late$max_day%in%myrange,])[1]
  resN=rbind(resN, data.frame(week.setting=ii,n.early=n.early, n.late=n.late))
}

####
junk=NULL
for(myday in days.range){
dat.lab.tmp=data_lab_clean2(code.dict, dat.x.raw, nm.value="value", day=myday)
dat.lab.tmp$AA=dat.lab.tmp$aspartate_aminotransferase_AST/dat.lab.tmp$alanine_aminotransferase_ALT

junk[[myday]]=dat.lab.tmp
}

dat.fit=NULL
for(nm.lab in nm.lab.keep){
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
tmp=left_join(tmp, dat.cc.stay, by="patient_num")
lab.max=data.table(tmp[complete.cases(tmp),c("day", "patient_num")])
lab.max=data.frame(lab.max[,max(day), by="patient_num"])
colnames(lab.max)[2]="max_lab_day"
tmp=left_join(tmp, lab.max, by="patient_num")
tmp=left_join(tmp,dat.cls, by="patient_num")
tmp$charlson_score=log(tmp$charlson_score+0.5)
dat.fit[[nm.lab]]=tmp
}

res=NULL
for(nm.lab in nm.lab.keep){
  tryCatch({
  print(nm.lab)
  dat.lab=dat.fit[[nm.lab]]
  dat.lab=dat.lab[complete.cases(dat.lab),]
  dat.lab=dat.lab[dat.lab$y%in%c(Inf, -Inf)!=1,]
  dat.lab=dat.lab[dat.lab$patient_num%in%patient_num.keep,]
  day.range.list=list()
  day.range.list[["week1"]]=c(0:7)
  day.range.list[["week2"]]=c(8:14)
  day.range.list[["week3"]]=c(15:60)
  day.range.list[["all"]]=c(0:60)
  
  for(ii in names(day.range.list)){
    myrange=day.range.list[[ii]]
    

 if(dim(dat.lab)[1]>10){
  if(is.ns==0){
  m1.early <- tryCatch(lmer(y ~ day+charlson_score+(day|patient_num), data=dat.lab[dat.lab$calendar_date<calendar.date.cut & dat.lab[,strat]%in%myrange,]),error=function(e) NA)
  m1.late <- tryCatch(lmer(y ~ day+charlson_score+(day|patient_num), data=dat.lab[dat.lab$calendar_date>=calendar.date.cut & dat.lab[,strat]%in%myrange,]),error=function(e) NA)

  s1.early=NA
  s1.late=NA
  ff=NA
  }else{
    
  tmp.range=unique(dat.lab[dat.lab[,strat]%in%myrange,"day"])
  tmp.range=intersect(tmp.range, days.range)
  if(length(tmp.range)<=7){myknots=myrange[c(3,5)]}
  if(length(tmp.range)>8){myknots=seq(min(tmp.range)+2, max(tmp.range)-2, by=3)}
    
  myknots=myknots[myknots<=max(dat.lab[dat.lab$calendar_date<calendar.date.cut& dat.lab[,strat]%in%myrange,"day"])]
  m1.early<- tryCatch(lmer(y ~ ns(day,knots=myknots)+charlson_score+(day|patient_num), data=dat.lab[dat.lab$calendar_date<calendar.date.cut& dat.lab[,strat]%in%myrange,]),error=function(e) NA)
  m1.late <- tryCatch(lmer(y ~ ns(day,knots=myknots)+charlson_score+(day|patient_num), data=dat.lab[dat.lab$calendar_date>=calendar.date.cut& dat.lab[,strat]%in%myrange,]),error=function(e) NA)

  ff=ns(c(days.range,dat.lab[dat.lab$calendar_date<calendar.date.cut& dat.lab[,strat]%in%myrange,]$day),knots=myknots)[1:length(days.range),]

  s1.early=tryCatch(cbind(1,ff, 0)%*%coef(summary(m1.early))[, "Estimate"],error=function(e) NA)
  s1.late=tryCatch(cbind(1,ff,0)%*%coef(summary(m1.late))[, "Estimate"], error=function(e) NA)
  }
 # plot: a1f1+a2f2+a3f3
 # a1,a2,a3 covariance matrix
 # knots: 7,14,21; cubic spline
 # day up to 30 days
 # risk prediction
 # 
  n.early=length(unique(dat.lab[dat.lab$calendar_date<calendar.date.cut,"patient_num"]))
  n.late=length(unique(dat.lab[dat.lab$calendar_date>=calendar.date.cut,"patient_num"]))
  
  n.lab.early=length(unique(dat.lab[dat.lab$calendar_date<calendar.date.cut& dat.lab[,strat]%in%myrange,"patient_num"]))
  n.lab.late=length(unique(dat.lab[dat.lab$calendar_date>=calendar.date.cut& dat.lab[,strat]%in%myrange,"patient_num"]))
  
  
  ss.early=tryCatch(sd(ranef(m1.early)$patient_num[,2]),error=function(e) NA)
  mm.early=tryCatch(coef(summary(m1.early))[-1, "Estimate"],error=function(e) NA)
  vcov.early=tryCatch((summary(m1.early))$vcov, error=function(e) NA)

  ss.late=tryCatch(sd(ranef(m1.late)$patient_num[,2]),error=function(e) NA)
  mm.late=tryCatch(coef(summary(m1.late))[-1, "Estimate"],error=function(e) NA)
  vcov.late=tryCatch((summary(m1.late))$vcov, error=function(e) NA)
  
  randomeffect.early=tryCatch(hist(ranef(m1.early)$patient_num[,2], breaks=100, freq=0), error=function(e) NA)
  randomeffect.late=tryCatch(hist(ranef(m1.late)$patient_num[,2], breaks=100, freq=0),error=function(e) NA)
  
  res[[nm.lab]][[ii]]=list(
    n.early=n.early,
    n.late=n.late,
    n.lab.early=n.lab.early,
    n.lab.late=n.lab.late,
    ss.early=ss.early, 
    mm.early=mm.early, 
    vcov.early=vcov.early,
    s1.early=s1.early,
    ss.late=ss.late, 
    mm.late=mm.late, 
    vcov.late=vcov.late,
    s1.late=s1.late,
    ff=ff
    #randomeffect.early=randomeffect.early, 
    #randomeffect.late=randomeffect.late
    )
 }
  }
  },error=function(e) NA)
}
res$resN=resN
res$summary.stay=summary.stay

res
}
