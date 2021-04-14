devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE, ref="master", force=T)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)
runAnalysis()

#setwd(.../FourCePhase2.1Survival/)
currSiteId = FourCePhase2.1Data::getSiteId()
load("data/betahat.port.rda")
load("data/code.dict.rda")

files.sources = list.files("R/utilis")
lapply(files.sources, function(xx) source(xx))

