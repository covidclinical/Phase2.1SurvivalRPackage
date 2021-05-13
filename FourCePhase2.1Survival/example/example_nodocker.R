rm(list=ls())
remove.packages("FourCePhase2.1Survival")
devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE, ref="master", force=T)
currSiteId = "MGB" ## change to your siteid
dir.input="/Users/chuanhong/Documents/Input"
dir.output="/Users/chuanhong/Documents/Output"
FourCePhase2.1Survival::runAnalysis_nodocker(currSiteId, dir.input, dir.output)
