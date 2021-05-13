betahat.port.deceased=NULL
site.list=list.files("../../Phase2.1SurvivalRSummariesPublic/")
site.list=site.list[grepl("_Result_R1.Rdata", site.list)]
site.list=gsub("_Result_R1.Rdata","",site.list)

for(mymodel in c("Lit3.DemCls", "LabCommon.DemCls")){
  for(submodel in c("impute")){
for(siteid in site.list){
  load(paste0("../../Phase2.1SurvivalRSummariesPublic/", siteid,"_Result_R1.Rdata"))
  betahat.port.deceased[[mymodel]][[submodel]][[siteid]]=get(paste0("survfit.coxnet.", mymodel, ".", submodel))$deceased$betahat
}
}
}
save(betahat.port.deceased, file="data/betahat.port.deceased.rda")