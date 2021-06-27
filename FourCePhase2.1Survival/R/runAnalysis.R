
runAnalysis=function(){
devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)
library(mice)
currSiteId = getSiteId()
tryCatch(FourCePhase2.1Data::runQC(currSiteId),error=function(e) print(e))
sink(file = file.path(getProjectOutputDirectory(), paste0(currSiteId, "_log.txt")), split = TRUE, append = FALSE)
cat("TEMPORAL TREND ANALYSIS\n")
tryCatch(runAnalysis_TemporalTrend(currSiteId),error=function(e) print(e))
cat("Phase1.1 PAPER REVISION\n")
tryCatch(runAnalysis_R1(currSiteId),error=function(e) print(e))
cat("DONE")
sink(file=NULL)
}


