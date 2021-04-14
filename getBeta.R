betahat.port=NULL
for(mymodel in c("Lit3.DemCls", "LabCommon.DemCls")){
  for(submodel in c("impute", "ind")){
for(siteid in c("FRBDX","MGB", "BIDMC", "NWU", "upenn", paste0("VA", c(1:5)))){
  load(paste0("../Phase2.1SurvivalRSummariesPublic/Input/res-R1-", siteid,".Rdata"))
  betahat.port[[mymodel]][[submodel]][[siteid]]=get(paste0("survfit.coxnet.", mymodel, ".", submodel))$betahat
}
}
}
save(betahat.port, file="FourCePhase2.1Survival/data/betahat.port.rda")