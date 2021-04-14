rm(list=ls())
devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE, ref="testing", force=T)
runAnalysis()
submitAnalysis()