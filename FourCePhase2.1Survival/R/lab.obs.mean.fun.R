lab.obs.mean.fun=function(code.dict, dat.x.raw, dat.calendar, nm.patient_num, nm.days_since_admission){
  myday.list=list(0, 1, c(0:1), c(0:7), c(0:max(dat.x.raw$days_since_admission, na.rm=T)))
  res.obs.mean.all=NULL
  for(ii in 1:length(myday.list)){
    myday=myday.list[[ii]]
    dat.lab.tmp=data_lab_clean2(code.dict, dat.x.raw, nm.value="value", day=myday) 
    dat.lab.tmp=left_join(dat.calendar, dat.lab.tmp, by="patient_num")                          
    cmonth.list=sort(unique(substr(dat.lab.tmp$calendar_date,1,7)))
    
    res.obs.mean=NULL
    count=0
    for(cmonth in cmonth.list){
      count=count+1
      dat.lab.tmp.sub=dat.lab.tmp[substr(dat.lab.tmp$calendar_date,1,7)==cmonth, ]
      res.obs.mean[[count]]=c(cmonth,dim(dat.lab.tmp.sub)[1], apply(dat.lab.tmp.sub[,-c(1:3)],2, function(x) mean(x, na.rm=T)))
    }
    res.obs.mean=do.call(rbind, res.obs.mean)
    res.obs.mean=data.frame(res.obs.mean)
    res.obs.mean[,-1]=apply(res.obs.mean[,-1],2, function(xx) as.numeric(as.character(xx)))
    colnames(res.obs.mean)[1:2]=c("calendar_date","n")
    res.obs.mean.all[[ii]]=res.obs.mean
  }
  res.obs.mean.all
}
