
runAnalysis_nodocker=function(siteid, dir.input, dir.output){
currSiteid=siteid
devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)
library(mice)
tryCatch(FourCePhase2.1Data::runQC_nodocker(currSiteId, dir.input),error=function(e) print(e))
sink(file = file.path(dir.output, paste0(currSiteId, "_log.txt")), split = TRUE, append = FALSE)
#cat("TEMPORAL TREND ANALYSIS\n")
#tryCatch(runAnalysis_TemporalTrend_nodocker(currSiteId, dir.input, dir.output),error=function(e) print(e))
#cat("Phase1.1 PAPER REVISION\n")
#tryCatch(runAnalysis_R1_nodocker(currSiteId, dir.input, dir.output),error=function(e) print(e))
cat("maxmin \n")
tryCatch(runAnalysis_maxmin_nodocker(currSiteId, dir.input, dir.output),error=function(e) print(e))

#cat("aki \n")
#tryCatch(runAnalysis_aki_nodocker(currSiteId, dir.input, dir.output),error=function(e) print(e))

cat("distributed cox \n")
tryCatch(runAnalysis_dc_nodocker(currSiteId, dir.input, dir.output),error=function(e) print(e))

cat("distributed cox \n")
tryCatch(runAnalysis_TemporalTrend_additional_nodocker(currSiteId, dir.input, dir.output),error=function(e) print(e))

cat("DONE")
sink(file=NULL)
}


