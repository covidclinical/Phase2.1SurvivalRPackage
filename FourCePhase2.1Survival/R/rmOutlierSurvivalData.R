rmOutlierSurvivalData=function (dat.survival,patient_num.out) 
{
  dat.calendar=dat.survival$dat.calendar
  dat.analysis.severe=dat.survival$dat.analysis.severe
  dat.analysis.deceased=dat.survival$dat.analysis.deceased
  dat.analysis.severedeceased=dat.survival$dat.analysis.severedeceased
  dat.calendar=dat.calendar[which(dat.calendar$patient_num%in%patient_num.out!=1),]
  dat.analysis.severe=dat.analysis.severe[which(dat.analysis.severe$patient_num%in%patient_num.out!=1),]
  dat.analysis.deceased=dat.analysis.deceased[which(dat.analysis.deceased$patient_num%in%patient_num.out!=1),]
  dat.analysis.severedeceased=dat.analysis.severedeceased[which(dat.analysis.severedeceased$patient_num%in%patient_num.out!=1),]
  dat.survival.rmout=dat.survival
  dat.survival.rmout$dat.calendar=dat.calendar
  dat.survival.rmout$dat.analysis.severe=dat.analysis.severe
  dat.survival.rmout$dat.analysis.deceased=dat.analysis.deceased
  dat.survival.rmout$dat.analysis.severedeceased=dat.analysis.severedeceased
  dat.survival.rmout
}
