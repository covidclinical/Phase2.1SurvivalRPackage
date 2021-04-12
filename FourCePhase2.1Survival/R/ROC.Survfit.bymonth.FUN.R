ROC.Survfit.bymonth.FUN=function(dat.label, score, t0.all, nm.event, ipw, is.combine=F, is.sep=F){
  X=dat.label[,"days_since_admission"]
  D=dat.label[,nm.event]
  month.all=sort(unique(substr(dat.label$calendar_date,1,7)))
  if(is.combine==T){month.list=c(month.all[month.all<"2020-07"],"Since July")}else{month.list=month.all}

res=vector("list", length = length(t0.all))
names(res)=t0.all
for(tt in t0.all){
  Y=  I(X <= as.numeric(tt))*D
  if(ipw==T){Gt=WGT.CEN(X, D, tt)}else{Gt.t=NULL}

  for (imonth in month.list){
    id.month=which(grepl(imonth,dat.label$calendar_date))
    if(imonth=="Since July" & is.combine==T){
      id.month=which(substr(dat.label$calendar_date,1,7)>="2020-07")}

  if(is.sep==F){junk.roc=ROC.Est.FUN(Y[id.month],score[id.month,colnames(score)==tt],yy0=0.5,fpr0=seq(0,1,0.01))}else{
    junk.roc=ROC.Est.FUN(Y[id.month],score[id.month,2],yy0=0.5,fpr0=seq(0,1,0.01))
  }
  auc=junk.roc[1]
  roc=matrix(junk.roc[-1], ncol=6)[-1,]
  colnames(roc)=c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
  res[[as.character(tt)]][[imonth]]=list(auc=auc, roc=roc)
}
}
res
}

