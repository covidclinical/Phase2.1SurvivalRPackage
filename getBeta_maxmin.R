setwd("/Users/chuanhong/Documents/GitHub/Phase2.1SurvivalRPackage/FourCePhase2.1Survival/")
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
colnames(B)=ls(betahat.port.deceased$LabCommon.DemCls$impute)
B[is.na(B)]=0
B=B[,setdiff(colnames(B), "ICSM")]

site.list=list.files("/Users/chuanhong/Documents/GitHub/Phase2.1SurvivalRSummariesPublic/")
site.list=site.list[grepl("_Result_R1.Rdata", site.list)]
site.list=unique(gsub("_Result_R1.Rdata", "", site.list))

beta.maxmin=NULL
for(siteid in site.list){
  load(paste0("/Users/chuanhong/Documents/GitHub/Phase2.1SurvivalRSummariesPublic/",siteid,"_Result_R1.Rdata"))
  xx=survfit.coxnet.LabCommon.DemCls.impute$deceased$xx
  n=survfit.coxnet.LabCommon.DemCls.impute$deceased$N.train
  Sigma =xx/n
  nm.maxmin=intersect(colnames(Sigma), rownames(B))
  output= beta_star(B[nm.maxmin,], Sigma[nm.maxmin, nm.maxmin], delta=0.5)
  #### estimation of maximin
  beta<-output$beta.est
  beta.maxmin[[siteid]]=beta
}
save(beta.maxmin, file="data/beta.maxmin.rda")
