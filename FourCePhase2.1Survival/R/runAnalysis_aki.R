
runAnalysis_aki=function(currSiteId){
dir.input=FourCePhase2.1Data::getInputDataDirectoryName()
dir.output=getProjectOutputDirectory()
obfuscation.level=FourCePhase2.1Data::getObfuscation(toupper(currSiteId))
obfuscation=F
if(length(obfuscation.level)!=0){
  if(obfuscation.level!=0){obfuscation=T}}
#### 1. read data
cat("1. read data\n")
LocalPatientObservations=FourCePhase2.1Data::getLocalPatientObservations_nodocker(currSiteId, dir.input)
LocalPatientClinicalCourse=FourCePhase2.1Data::getLocalPatientClinicalCourse_nodocker(currSiteId, dir.input)
LocalPatientSummary=FourCePhase2.1Data::getLocalPatientSummary_nodocker(currSiteId, dir.input)
data(code.dict, package="FourCePhase2.1Data")
cat("data pivot \n")


month.list=list(
paste0("2020-0", c(1:4)),
paste0("2020-0", c(5:6)),
paste0("2020-0", c(7:8)),
c("2020-09","2020-10"),
c(paste0("2020-", c(11:12)),"2021-01"),
paste0("2020-0", c(1:6)),
c(paste0("2020-0", c(7:9)), paste0("2020-", c(10:12)),"2021-01")
)

res=vector(mode = "list", length = length(month.list))
res.nm=unlist(lapply(month.list, function(xx){xx1=xx[1]; xx2=xx[length(xx)];paste(xx1, xx2, sep=":")}))
names(res)=res.nm
for(ll in 1:length(month.list)){
demog=LocalPatientSummary[substr(LocalPatientSummary$admission_date,1,7)%in%month.list[[ll]],]
obs=LocalPatientObservations[LocalPatientObservations$patient_num%in%demog$patient_num,]
tmp=get_aki_events(demog,obs)
res[[res.nm[ll]]]=tmp
}
res.aki=res
save(res.aki=res.aki,
     file=file.path(dir.output, paste0(currSiteId, "_aki.Rdata")))
}


