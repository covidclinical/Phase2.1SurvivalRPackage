version_control_fun=function(dir.input.old="/Users/chuanhong/Documents/Input/", dir.input.new="/Users/chuanhong/Documents/Input_new/"){
dat.survival=getSurvivalData(dir.input.old, code.dict, siteid="MGB")
dat.survival.new=getSurvivalData(dir.input.new, code.dict, siteid="MGB")

table(dat.survival$dat.analysis.severe$severe)
table(dat.survival.new$dat.analysis.severe$severe)

table(dat.survival$dat.analysis.severedeceased$severedeceased)
table(dat.survival.new$dat.analysis.severedeceased$severedeceased)

table(dat.survival$dat.analysis.deceased$deceased)
table(dat.survival.new$dat.analysis.deceased$deceased)

xx0=read.csv(paste0(dir.input.old,"/LocalPatientClinicalCourse.csv"))
xx=xx0[xx0$days_since_admission==30,c("patient_num", "deceased")]
table(xx[,2])

yy0=read.csv(paste0(dir.input.new, "/LocalPatientClinicalCourse.csv"))
yy=yy0[yy0$days_since_admission==30,c("patient_num", "deceased")]
table(yy[,2])
table(unique(xx0$patient_num)%in%unique(yy0$patient_num))
}