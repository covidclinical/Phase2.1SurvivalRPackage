# Phase2.1SurvivalRPackage
R code to run, validate, and submit the analysis for the Survival project.

To install this package in R:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)
currSiteId = getSiteId()
runAnalysis(currSiteId)
submitAnalysis()
```


