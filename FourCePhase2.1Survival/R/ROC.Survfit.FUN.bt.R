ROC.Survfit.FUN.bt=function(dat.label, score, t0.all, nm.event, ipw, is.sep){
  dat.label=left_join(data.frame(patient_num=score$patient_num), dat.label, by="patient_num")
  X=dat.label[,"days_since_admission"]
  D=dat.label[,nm.event]
res=vector("list", length = length(t0.all))
names(res)=t0.all
names(res)=t0.all
for(tt in t0.all){
  Y=  I(X <= as.numeric(tt))*D
  if(ipw==T){Gt=WGT.CEN(X, D, as.numeric(tt))}else{Gt.t=NULL}
  if(is.sep==F){
  junk.roc=ROC.Est.FUN(Y,score[,colnames(score)==tt],yy0=0.5,fpr0=seq(0,1,0.01))}else{
  junk.roc=ROC.Est.FUN(Y,score[,2],yy0=0.5,fpr0=seq(0,1,0.01))
  }
  
  
  auc=junk.roc[1]
  roc=matrix(junk.roc[-1], ncol=6)[-1,]
  colnames(roc)=c("cut", "p.pos", "fpr", "tpr", "ppv", "npv")
  res[[as.character(tt)]]=list(auc=auc, roc=roc)
}
res
}

