lab.obs.mean.std.fun=function(code.dict, dat.x.raw, dat.calendar, nm.patient_num, nm.days_since_admission){
  myday.list=list(0, 1, c(0:1), c(0:7), c(0:max(dat.x.raw$days_since_admission, na.rm=T)))
  res.obs.mean.std.all=NULL
  for(ii in 1:length(myday.list)){
    myday=myday.list[[ii]]
    dat.lab.tmp=data_lab_clean2(code.dict, dat.x.raw, nm.value="value", day=myday) 
    dat.lab.tmp=left_join(dat.calendar, dat.lab.tmp, by="patient_num")                          
    cmonth.list=sort(unique(substr(dat.lab.tmp$calendar_date,1,7)))
    
    res.obs.mean.std=NULL
    count=0
    sd.all=apply(dat.lab.tmp[,-c(1:3)],2, function(x) sd(x, na.rm=T))
    for(cmonth in cmonth.list){
      count=count+1
      dat.lab.tmp.sub=dat.lab.tmp[substr(dat.lab.tmp$calendar_date,1,7)==cmonth, ]
      junk=apply(dat.lab.tmp.sub[,-c(1:3)],2, function(x) mean(x, na.rm=T))/sd.all
      res.obs.mean.std[[count]]=c(cmonth,dim(dat.lab.tmp.sub)[1],junk)
    }
    res.obs.mean.std=do.call(rbind, res.obs.mean.std)
    res.obs.mean.std=data.frame(res.obs.mean.std)
    res.obs.mean.std[,-1]=apply(res.obs.mean.std[,-1],2, function(xx) as.numeric(as.character(xx)))
    colnames(res.obs.mean.std)[1:2]=c("calendar_date","n")
    res.obs.mean.std.all[[ii]]=res.obs.mean.std
  }
  res.obs.mean.std.all
}
