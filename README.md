# Phase2.1SurvivalRPackage for Docker user.
R code to run, validate, and submit the analysis for the Phase2.1 Temporal Trend project and additional results for Phase1.1 paper revision. Please note it will take about 30 minutes to finish running the analysis. Thank you for your patience!

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

Share with @Chuan your GitHub handle via direct message on Slack channel so you can be added as contributor to the repository. Note that you would need to use a token to access private repos, see here. https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token

Briefly, to generate a new token, go to your GitHub settings -> Developer settings -> Personal access tokens -> Generate.

If somehow submitAnalysis() didn’t allow you to upload the results to Phase2.1SurvivalRSummariesPublic, you can share the results file with @Chuan via the direct message on Slack channel.
