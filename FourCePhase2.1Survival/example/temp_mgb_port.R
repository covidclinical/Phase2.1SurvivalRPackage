
survfit.coxnet.port.betahat=NULL
data(betahat.port.deceased, package="FourCePhase2.1Survival")
submodel="impute"
for(mysite in ls(betahat.port.deceased$Lit3.DemCls$impute)){
  print(mysite)
  for(mymodel in c("Lit3.DemCls", "LabCommon.DemCls")){
    betahat=betahat.port.deceased[[mymodel]][[submodel]][[mysite]]
    survfit.coxnet.port.betahat[[mymodel]][[submodel]][[mysite]]=tryCatch(survfit.glmnet.coefficient.R1.fun(dat.survival, ipw=T, nm.event="deceased", nm.lab.all=nm.lab.LabAll, betahat= betahat, nm.cls, siteid, dir.output, 
                                                                                                            period.train, period.valid, calendar.date.cut="2020-07",  t0.all=c(1:14), yes.cv=F, is.bt=T),error=function(e){print(e); NA})
  }
}
survfit.coxnet.port.betahat.mgb=survfit.coxnet.port.betahat
save(survfit.coxnet.port.betahat.mgb, file="/4ceData/GitHub/Phase2.1SurvivalRSummariesPublic/Output_R1/MGB_port_temp.Rdata")