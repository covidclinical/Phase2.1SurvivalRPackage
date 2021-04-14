rm(list=ls())
devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE, ref="master", force=T)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)

currSiteId = "MGB" ## change to your siteid
dir.input="/Users/chuanhong/Documents/Input"
dir.output="/Users/chuanhong/Documents/Output"

runAnalysis_nodocker(currSiteId, dir.input, dir.output)
