runAnalysis_VA=function(currSiteId, dir.output){
  file.log=file(file.path(dir.output, paste0(currSiteId, "_log.txt")),open = "wt")
  sink(file.log)
  sink(file.log, type="message")
  cat("TEMPORAL TREND ANALYSIS\n")
  tryCatch(runAnalysis_TemporalTrend_VA(currSiteId, dir.package, dir.output),error=function(e) print(e))
  cat("Phase1.1 PAPER REVISION\n")
  tryCatch(runAnalysis_R1_VA(currSiteId, dir.package, dir.output),error=function(e) print(e))
  sink(type = "message")
  sink()
}