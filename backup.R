
#### 11. Federated learning
cat("11.federated learning")
survfit.coxnet.port.betahat.deceased=NULL
data(betahat.port.deceased, package="FourCePhase2.1Survival")
submodel="impute"

B.nm=unique(unlist(lapply(ls(betahat.port.deceased$LabCommon.DemCls$impute), function(ll) names(betahat.port.deceased[["LabCommon.DemCls"]][[submodel]][[ll]]))))

B=lapply(ls(betahat.port.deceased$LabCommon.DemCls$impute), function(ll){
  B.junk=betahat.port.deceased[["LabCommon.DemCls"]][[submodel]][[ll]]
  B.junk=left_join(data.frame(var=B.nm), data.frame(var=names(B.junk),coefficient=B.junk), by="var")
  B.junk=B.junk$coefficient
  B.junk
}
)
B=do.call(cbind, B)
rownames(B)=B.nm
B[is.na(B)]=0

junk.xx=getX.coxnet.R1.fun(dat.survival, nm.event="deceased", nm.lab.keep=nm.9lab, nm.dem, nm.cls, siteid, dir.output, 
                           period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=T, K=10, is.bt=T, include.ind=F)$dat
mydat=junk.xx[,B.nm]
xx=t(data.matrix(mydat))%*%data.matrix(mydat)

Sigma =xx /dim(mydat)[1]
output= beta_star(B, Sigma, delta=0.5)
#### estimation of maximin
beta<-output$beta.est
#### maximin weight
output$weight
####### Check Accuracy on each group ########
accuracy <- function(X, Y, b){
  logit = exp(X%*%b) / (1+exp(X%*%b))
  Y_hat = as.integer(logit >= 0.5)
  return(mean(Y_hat==Y)) 
}

for(mysite in ls(betahat.port.deceased$LabCommon.DemCls$impute)){
  print(mysite)
  survfit.coxnet.port.betahat.deceased[[mymodel]][[submodel]][[mysite]]=tryCatch(survfit.glmnet.coefficient.R1.fun(dat.survival, ipw=T, nm.event="deceased", nm.lab.all=nm.lab.LabAll, betahat= betahat, nm.cls, siteid, dir.output, 
                                                                                                                   period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=F, is.bt=T),error=function(e){print(e); NA})
}
