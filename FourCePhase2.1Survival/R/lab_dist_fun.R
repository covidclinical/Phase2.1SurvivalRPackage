lab_dist_fun=function(dir.input, code.dict, LocalPatientObservations, dat.survival,calendar.date.cut="2020-07", myscale, lab.breaks.original,lab.breaks.log){
dat.calendar=dat.survival$dat.calendar
dat.x.raw =LocalPatientObservations
patient_num.severe=dat.survival$dat.analysis.severe$patient_num[which(dat.survival$dat.analysis.severe$severe==1)]
patient_num.deceased=dat.survival$dat.analysis.deceased$patient_num[which(dat.survival$dat.analysis.deceased$deceased==1)]
patient_num.severedeceased=dat.survival$dat.analysis.severedeceased$patient_num[which(dat.survival$dat.analysis.severedeceased$severedeceased==1)]

res.all=obs.all=NULL
for(nm.event in c("all", "deceased", "severedeceased")){
res=obs=NULL  
for(myday.list in c(1:4)){
  if(myday.list==1){myday=0}
  if(myday.list==2){myday=1}
  if(myday.list==3){myday=c(0,1)}
  if(myday.list==4){myday=c(0:7)}
  
if(nm.event=="all"){
tmp=data_lab_clean2(code.dict, dat.x.raw, nm.value="value", day=myday) 
}
  
if(nm.event=="severedeceased"){
tmp=data_lab_clean2(code.dict, dat.x.raw[dat.x.raw$patient_num%in%patient_num.severedeceased,], nm.value="value", day=myday) 
}
  
if(nm.event=="deceased"){
tmp=data_lab_clean2(code.dict, dat.x.raw[dat.x.raw$patient_num%in%patient_num.deceased,], nm.value="value", day=myday) 
}
tmp$AA=tmp$aspartate_aminotransferase_AST/tmp$alanine_aminotransferase_ALT
tmp=left_join(tmp, dat.calendar,by="patient_num")
nm.lab.all=setdiff(colnames(tmp), c("patient_num", "calendar_date", "calendar_day"))
lab.breaks=lab.breaks.original
lab.breaks=rbind(lab.breaks, data.frame(labname="AA", break.lb=0, break.ub=max(tmp$AA,na.rm=T)))
if(myscale=="log"){
tmp[,nm.lab.all]=log(tmp[,nm.lab.all])
lab.breaks=lab.breaks.log
}
tmp.early=tmp[tmp$calendar_date<calendar.date.cut,]
tmp.late=tmp[tmp$calendar_date>=calendar.date.cut,]

res.early=lapply(nm.lab.all, function(xx) tryCatch(hist(tmp.early[,xx],
                                                        breaks=c(lab.breaks$break.lb[lab.breaks$labname==xx]-10, seq(lab.breaks$break.lb[lab.breaks$labname==xx], lab.breaks$break.ub[lab.breaks$labname==xx],(lab.breaks$break.ub[lab.breaks$labname==xx]-lab.breaks$break.lb[lab.breaks$labname==xx])/100),lab.breaks$break.ub[lab.breaks$labname==xx]*100), 
                                                        na.rm=T, plot=F),error=function(e) NA))
res.late=lapply(nm.lab.all, function(xx) tryCatch(hist(tmp.late[,xx],
                                                       breaks=c(lab.breaks$break.lb[lab.breaks$labname==xx]-10, seq(lab.breaks$break.lb[lab.breaks$labname==xx], lab.breaks$break.ub[lab.breaks$labname==xx],(lab.breaks$break.ub[lab.breaks$labname==xx]-lab.breaks$break.lb[lab.breaks$labname==xx])/100),lab.breaks$break.ub[lab.breaks$labname==xx]*100), 
                                                       na.rm=T, plot=F), error=function(e) NA))
res.early=lapply(nm.lab.all, function(xx) tryCatch(hist(tmp.early[,xx],
                                                        breaks=c(lab.breaks$break.lb[lab.breaks$labname==xx]-10, seq(lab.breaks$break.lb[lab.breaks$labname==xx], lab.breaks$break.ub[lab.breaks$labname==xx],(lab.breaks$break.ub[lab.breaks$labname==xx]-lab.breaks$break.lb[lab.breaks$labname==xx])/100),lab.breaks$break.ub[lab.breaks$labname==xx]*100), 
                                                        na.rm=T, plot=F),error=function(e) NA))
res.both=lapply(nm.lab.all, function(xx) tryCatch(hist(tmp[,xx],
                                                       breaks=c(lab.breaks$break.lb[lab.breaks$labname==xx]-10, seq(lab.breaks$break.lb[lab.breaks$labname==xx], lab.breaks$break.ub[lab.breaks$labname==xx],(lab.breaks$break.ub[lab.breaks$labname==xx]-lab.breaks$break.lb[lab.breaks$labname==xx])/100),lab.breaks$break.ub[lab.breaks$labname==xx]*100), 
                                                       na.rm=T, plot=F), error=function(e) NA))

names(res.early)=names(res.late)=names(res.both)=nm.lab.all
res[[myday.list]]=list(res.early=res.early,res.late=res.late, res.both=res.both)

obs.early=lapply(nm.lab.all, function(xx) sum(is.na(tmp.early[,xx])!=1))
obs.late=lapply(nm.lab.all, function(xx) sum(is.na(tmp.late[,xx])!=1))
obs.both=lapply(nm.lab.all, function(xx) sum(is.na(tmp[,xx])!=1))

names(obs.early)=names(obs.late)=names(obs.both)=nm.lab.all
obs[[myday.list]]=list(obs.early=obs.early,obs.late=obs.late, obs.both=obs.both)
}
names(res)=names(obs)=c(1:4)
res.all[[nm.event]]=res
obs.all[[nm.event]]=obs
}
list(res.all=res.all, obs.all=obs.all)
}
