
runAnalysis=function(){
file.log=file(file.path(getProjectOutputDirectory(), paste0(currSiteId, "_log.txt")),open = "wt")
sink(file.log)
sink(file.log, type="message")
cat("TEMPORAL TREND ANALYSIS\n")
tryCatch(runAnalysis_TemporalTrend(currSiteId),error=function(e) print(e))
cat("Phase1.1 PAPER REVISION\n")
tryCatch(runAnalysis_R1(currSiteId),error=function(e) print(e))
sink(type = "message")
sink()
}


