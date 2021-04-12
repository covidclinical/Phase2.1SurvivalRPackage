baseline.score.lab.fun=function(dat.survival, siteid, dir.output){
  mycol=c("red", "blue", "green")
  mycol.lt <- c(t_col("red", perc = 70, name = "lg.red"),
                t_col("blue", perc = 70, name = "lg.blue"),
                t_col("green", perc = 70, name = "lg.green"))
  
  nm.event="severedeceased"
  dat=dat.survival[[paste0("dat.analysis.",nm.event)]]
  dat=dat[,duplicated(colnames(dat))!=1]
  dat.calendar=dat.survival[["dat.calendar"]]
  dat=left_join(dat.calendar[,c("patient_num", "calendar_day", "calendar_date")], dat, by="patient_num")
  dim(dat)
  dat=dat[,c(nm.patient_num, nm.days_since_admission, nm.event, "calendar_date","calendar_day", nm.dem, nm.lab.keep)]
  dat$age_group_new = relevel(dat$age_group_new, ref = "50to69")
  dat$sex=tolower(dat$sex)
  dat$sex=factor(dat$sex, level=c("male", "female"))
  dat$race_new=factor(dat$race_new, level=c("White","Black","Asian","Hispanic and Other"))
  dat.lab.impute=do.call(cbind,lapply(nm.lab.keep, function(ll){x=dat[,ll];mm=mean(x,na.rm=T);x[is.na(x)]=mm; x}))
  obs.indx=1*do.call(cbind,lapply(nm.lab.keep, function(ll){x=dat[,ll];is.na(x)!=1}))
  colnames(obs.indx)=paste0("obs_", nm.lab.keep)
  dat[,nm.lab.keep]=dat.lab.impute
  dat=cbind(dat, obs.indx)
  colnames(dat)=gsub("C_reactive_protein_CRP_Normal_Sensitivity", "CRP", colnames(dat))
  colnames(dat)=gsub("lymphocyte_count", "LYM", colnames(dat))
  colnames(dat)=gsub("D_dimer", "DD", colnames(dat))
  colnames(dat)=gsub("age_group_new", "age", colnames(dat))
  colnames(dat)=gsub("race_new", "race", colnames(dat))
  head(dat)
  dat$DD=log(dat$DD+1)
  
  ##### 1. risk score and lab measures at day 0
  dat.day0=dat
  Y=ifelse(dat$days_since_admission==0 & dat[,nm.event]==1 ,1,0)
  
  multi.formulas = as.formula(paste('Y~',
                                    paste(nm.dem.new, collapse="+"),"+",
                                    paste(nm.lab.new,collapse="+"), "+",
                                    paste(paste0("obs_", nm.lab.new),collapse="+")))
  
  dat.day0.new=model.matrix(multi.formulas,dat.day0)[,-1]
  junk=Est.ALASSO.GLM(cbind(Y, dat.day0.new), regularize=T)
  score.day0=g.logit(cbind(1,dat.day0.new)%*%junk[1:(dim(dat.day0.new)[2]+1)])
  junk.day0=ROC.Est.FUN(Y,score.day0,yy0=0.5,fpr0=seq(0,1,0.01),wgti=NULL,yes.smooth=F)
  mtx.day0=matrix(junk.day0[-1], ncol=6)
  colnames(mtx.day0)=c("cut","p.pos","fpr","tpr","ppv","npv")
  junk1=mtx.day0[which.min(abs(mtx.day0[,"fpr"]-0.1)),]
  junk2=mtx.day0[which.min(abs(mtx.day0[,"tpr"]-0.9)),]
  
  cut1=junk1["cut"]
  cut2=junk2["cut"]
  
  score.bin.day0=rep("M",dim(dat)[1])
  score.bin.day0[which(score.day0<cut2)]="L"
  score.bin.day0[which(score.day0>cut1)]="H"
  
  ppv.H=mean(dat.day0$severedeceased[score.bin.day0=="H"])
  ppv.M=mean(dat.day0$severedeceased[score.bin.day0=="M"])
  ppv.L=mean(dat.day0$severedeceased[score.bin.day0=="L"])
  
  cday.list=sort(unique(dat.day0$calendar_day))
  junk.score=lapply(cday.list, function(ll){
    temp=rep(0,3); names(temp)=c("H", "M", "L")
    res.temp=prop.table(table(score.bin.day0[dat.day0$calendar_day==ll]))
    temp[names(res.temp)]=res.temp
    c(cday=ll, temp)})
  res.score=data.frame(do.call(rbind, junk.score))
  
  cmonth.list=sort(unique(substr(dat.day0$calendar_date,1,7)))
  junk.score.month=lapply(cmonth.list, function(ll){
    temp=rep(0,3); names(temp)=c("H", "M", "L")
    res.temp=prop.table(table(score.bin.day0[substr(dat.day0$calendar_date,1,7)==ll]))
    temp[names(res.temp)]=res.temp
    c(cmonth=ll, temp)})
  res.score.month=data.frame(do.call(rbind, junk.score.month))
  res.score.month.trans=data.frame(t(res.score.month[,-1]))
  colnames(res.score.month.trans)=res.score.month[,1]
  res.score.month.trans=apply(res.score.month.trans,2, as.numeric)
  
  
  
  junk.lab.month=lapply(cmonth.list, function(ll){
    temp=rep(0,3); names(temp)=c("obs_DD", "obs_CRP", "obs_LYM")
    temp=c(mean(dat.day0$obs_DD[substr(dat.day0$calendar_date,1,7)==ll]),mean(dat.day0$obs_CRP[substr(dat.day0$calendar_date,1,7)==ll]),mean(dat.day0$obs_LYM[substr(dat.day0$calendar_date,1,7)==ll]))
    c(cmonth=ll, temp)})
  res.lab.month=data.frame(do.call(rbind, junk.lab.month))
  res.lab.month.trans=data.frame(t(res.lab.month[,-1]))
  colnames(res.lab.month.trans)=res.lab.month[,1]
  res.lab.month.trans=apply(res.lab.month.trans,2, as.numeric)
  rownames(res.lab.month.trans)=c("obs_DD", "obs_CRP", "obs_LYM")
  
  pdf(file=paste0(dir.output,"/bar.month.",nm.event,".", siteid,".pdf"), height=6, width=12)
  par(mfrow=c(1,1))
  barplot(res.score.month.trans, beside = TRUE, col=mycol, ylim=c(0,0.7))
  legend("topleft", c("H", "M", "L"), col=mycol, pch=15, ncol=3)
  
  #barplot(res.lab.month.trans, beside = TRUE, col=mycol, ylim=c(0,0.7))
  #legend("topleft", c("obs_DD", "obs_CRP", "obs_LYM"), col=mycol, pch=15, ncol=3)
  
  dev.off()
  
  junk.lab=lapply(cday.list, function(ll){
    temp.DD=mean(dat.day0$obs_DD[dat.day0$calendar_day==ll])
    temp.CRP=mean(dat.day0$obs_CRP[dat.day0$calendar_day==ll])
    temp.LYM=mean(dat.day0$obs_LYM[dat.day0$calendar_day==ll])
    c(cday=ll, obs.DD=temp.DD, obs.CRP=temp.CRP, obs.LYM=temp.LYM)})
  res.lab=data.frame(do.call(rbind, junk.lab))
  
  bw.score.H <- np::npregbw(formula =res.score$H~ res.score$cday, regtype="ll")
  bw.score.M <- np::npregbw(formula =res.score$M~ res.score$cday, regtype="ll")
  bw.score.L <- np::npregbw(formula =res.score$L~ res.score$cday, regtype="ll")
  
  kre.score.H <- np::npreg(bws = bw.score.H)
  kre.score.M <- np::npreg(bws = bw.score.M)
  kre.score.L <- np::npreg(bws = bw.score.L)
  
  bw.lab.DD <- np::npregbw(formula =res.lab$obs.DD~ res.lab$cday, regtype="ll")
  bw.lab.CRP <- np::npregbw(formula =res.lab$obs.CRP~ res.lab$cday, regtype="ll")
  bw.lab.LYM <- np::npregbw(formula =res.lab$obs.LYM~ res.lab$cday, regtype="ll")
  
  kre.lab.DD <- np::npreg(bws = bw.lab.DD)
  kre.lab.CRP <- np::npreg(bws = bw.lab.CRP)
  kre.lab.LYM <- np::npreg(bws = bw.lab.LYM)
  
  res.score$cdate=res.score$cday+as.Date("2020-03-11")
  res.lab$cdate=res.lab$cday+as.Date("2020-03-11")
  
  pdf(file=paste0(dir.output,"/score.calendar.",nm.event,".", siteid,".pdf"), height=6, width=12)
  par(mfrow=c(1,2))
  plot(res.score$cday, res.score$H, col=mycol.lt[1], xlab="", ylab="", main=paste0("Risk of being ", nm.event, " at day 0\n AUC=", round(junk.day0[1],2),"; ppv.H=", round(ppv.H,2), "; ppv.M=", round(ppv.M,2), "; ppv.L=", round(ppv.L,2)), xaxt="n", lty=2, pch=16)
  points(res.score$cday, res.score$M, col=mycol.lt[2],lty=2, pch=16)
  points(res.score$cday, res.score$L, col=mycol.lt[3], lty=2, pch=16)
  
  #lines(kre.score.H$eval$`res.score$cday`, kre.score.H$mean, type="l", col=mycol[1], lwd=2)
  #lines(kre.score.M$eval$`res.score$cday`, kre.score.M$mean, type="l", col=mycol[2], lwd=2)
  #lines(kre.score.L$eval$`res.score$cday`, kre.score.L$mean, type="l", col=mycol[3], lwd=2)
  
  res.score.ave=data.frame(do.call(rbind,lapply(6:max(res.score$cday), function(ll) c(cday=ll, colMeans(res.score[res.score$cday%in%c(ll:(ll+6)), c("H", "M", "L")], na.rm=T)))))
  res.score.ave$H[is.na(res.score.ave$H)]=mean(res.score.ave$H,na.rm=T)
  res.score.ave$M[is.na(res.score.ave$M)]=mean(res.score.ave$M,na.rm=T)
  res.score.ave$L[is.na(res.score.ave$L)]=mean(res.score.ave$L,na.rm=T)
  
  lines(res.score.ave$cday, res.score.ave$H, type="l", col=mycol[1], lwd=2)
  lines(res.score.ave$cday, res.score.ave$M, type="l", col=mycol[2], lwd=2)
  lines(res.score.ave$cday, res.score.ave$L, type="l", col=mycol[3], lwd=2)
  
  legend("topleft", c("H","M", "L"), lty=1, col=mycol, ncol=3, lwd=2)
  axis(side=1,at=which(as.character(res.score$cdate)%in%c("2020-03-11", "2020-04-01", "2020-05-01", "2020-06-01", "2020-07-15")), labels=F)
  text(x = which(as.character(res.score$cdate)%in%c("2020-03-11", "2020-04-01", "2020-05-01", "2020-06-01", "2020-07-15")),
       y = par("usr")[3] - 0.12,
       labels = c("2020-03-11","2020-04-01", "2020-05-01", "2020-06-01", "2020-07-15"),
       xpd = NA,
       srt = 35,
       cex = 1)
  
  plot(res.lab$cday, res.lab$obs.DD, col=mycol.lt[1], xlab="", xaxt="n", ylim=c(0,1),ylab="", main=paste0("Lab measurement rate at day 0"), lty=2, pch=16)
  points(res.lab$cday, res.lab$obs.CRP, col=mycol.lt[2], lty=2, pch=16)
  points(res.lab$cday, res.lab$obs.LYM, col=mycol.lt[3], lty=2, pch=16)
  
  #lines(kre.lab.DD$eval$`res.lab$cday`, kre.lab.DD$mean, type="l", col="red", lwd=2)
  #lines(kre.lab.CRP$eval$`res.lab$cday`, kre.lab.CRP$mean, type="l", col="blue", lwd=2)
  #lines(kre.lab.LYM$eval$`res.lab$cday`, kre.lab.LYM$mean, type="l", col="green", lwd=2)
  
  res.lab.ave=data.frame(do.call(rbind,lapply(6:max(res.lab$cday), function(ll) c(cday=ll, colMeans(res.lab[res.lab$cday%in%c(ll:(ll+6)), c("obs.DD", "obs.CRP", "obs.LYM")], na.rm=T)))))
  res.lab.ave$obs.DD[is.na(res.lab.ave$obs.DD)]=mean(res.lab.ave$obs.DD,na.rm=T)
  res.lab.ave$obs.CRP[is.na(res.lab.ave$obs.CRP)]=mean(res.lab.ave$obs.CRP,na.rm=T)
  res.lab.ave$obs.LYM[is.na(res.lab.ave$obs.LYM)]=mean(res.lab.ave$obs.LYM,na.rm=T)
  
  lines(res.lab.ave$cday, res.lab.ave$obs.DD, type="l", col=mycol[1], lwd=2)
  lines(res.lab.ave$cday, res.lab.ave$obs.CRP, type="l", col=mycol[2], lwd=2)
  lines(res.lab.ave$cday, res.lab.ave$obs.LYM, type="l", col=mycol[3], lwd=2)
  
  legend("topleft", c("DD","CRP", "LYM"), lty=1, col=c("red","blue","green"), ncol=3, lwd=2)
  #lines(kre.lab.DD$eval$`res.lab$cday`, kre.lab.DD$mean, type="l", col="red", lwd=2)
  #lines(kre.lab.CRP$eval$`res.lab$cday`, kre.lab.CRP$mean, type="l", col="blue", lwd=2)
  #lines(kre.lab.LYM$eval$`res.lab$cday`, kre.lab.LYM$mean, type="l", col="green", lwd=2)
  legend("topleft", c("DD","CRP", "LYM"), lty=1, col=c("red","blue","green"), ncol=3, lwd=2)
  axis(side=1,at=which(as.character(res.lab$cdate)%in%c("2020-03-11", "2020-04-01", "2020-05-01", "2020-06-01", "2020-07-15")), labels=F)
  text(x = which(as.character(res.lab$cdate)%in%c("2020-03-11", "2020-04-01", "2020-05-01", "2020-06-01", "2020-07-15")),
       y = par("usr")[3] - 0.12,
       labels = c("2020-03-11","2020-04-01", "2020-05-01", "2020-06-01", "2020-07-15"),
       xpd = NA,
       srt = 35,
       cex = 1)
  dev.off()
  dat.day0.save=dat.day0[,c("patient_num","days_since_admission", "severedeceased" , "calendar_date" , "calendar_day", "obs_DD", "obs_CRP", "obs_LYM")]
  dat.day0.save=data.frame(score=score.bin.day0, dat.day0.save)
  dat.day0.save
}