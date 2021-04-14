# Phase2.1SurvivalRPackage for Docker User, Non-Docker User, and VA User. 

# 1. Docker User
R code to run, validate, and submit the analysis for the Phase2.1 Temporal Trend project and additional results for Phase1.1 paper revision. Please note it will take about 30 minutes to finish running the analysis. Thank you for your patience!

## 1. Make sure you are using the latest version of Docker. 

## 2. Always RESTART your R session before installing or re-installing the package!

## 3. Run the following scripts in R:

```
devtools::install_github("https://github.com/covidclinical/Phase2.1SurvivalRPackage", subdir="FourCePhase2.1Survival", upgrade=FALSE)
library(FourCePhase2.1Data)
library(FourCePhase2.1Survival)
library(icd)
currSiteId = getSiteId()
runAnalysis(currSiteId)
```

## 4. Two ways to submit the results:
1. Submit via GitHub. 
+ (1) Share with @Chuan your GitHub accountname via direct message on Slack channel so you can be added as contributor to the repository. 
+ (2) Note that you would need to use a token to access private repos, see here (https://docs.github.com/en/github/authenticating-to-github/creating-a-personal-access-token). Briefly, to generate a new token, go to your GitHub settings -> Developer settings -> Personal access tokens -> Generate.
+ (3) Submit the result files in your output directory to the branch named topic-[YourSiteId] in GitHub repo Phase2.1SurvivalRSummariesPublic by running the following scripts in R:
```
submitAnalysis()
```

2. Submit via Slack channel. Alternatively, if somehow submitAnalysis() didn’t allow you to upload the results to Phase2.1SurvivalRSummariesPublic, you can share the results file with @Chuan via the direct message on Slack channel.
