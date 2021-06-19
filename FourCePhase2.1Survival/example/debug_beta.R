devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE, ref="master", force=T)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)
#runAnalysis()

#setwd("/4ceData/GitHub/Phase2.1SurvivalRpackage/FourCePhase2.1Survival/")
setwd("FourCePhase2.1Survival/")

#currSiteId = FourCePhase2.1Data::getSiteId()
currSiteId = "MGB" ## change to your siteid
dir.input="/Users/chuanhong/Documents/Input"
dir.output="/Users/chuanhong/Documents/Output"


load("data/betahat.port.deceased.rda")
load("data/betahat.port.rda")
load("data/code.dict.rda")

files.sources = list.files("R")
lapply(files.sources, function(xx) source(paste0("R/",xx)))

