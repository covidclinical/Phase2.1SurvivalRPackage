data_dem_clean=function (dat.dem.raw, nm.patient_num, nm.gender, nm.age, nm.race) 
{
  nm.age.num = "age_cat"
  dat.sub = data.frame(dat.dem.raw[, c(nm.patient_num, nm.gender, 
                                       nm.age, nm.race)])
  dat.sub[, nm.age] = as.character(dat.sub[, nm.age])
  dat.sub[, nm.gender] = as.character(dat.sub[, nm.gender])
  dat.sub[, nm.race] = as.character(dat.sub[, nm.race])
  dat.sub[, nm.race] = str_to_title(dat.sub[, nm.race])
  dat.sub[dat.sub[, nm.age] == "00to02", nm.age.num] = 1
  dat.sub[dat.sub[, nm.age] == "03to05", nm.age.num] = 2
  dat.sub[dat.sub[, nm.age] == "06to11", nm.age.num] = 3
  dat.sub[dat.sub[, nm.age] == "12to17", nm.age.num] = 4
  dat.sub[dat.sub[, nm.age] == "18to25", nm.age.num] = 5
  dat.sub[dat.sub[, nm.age] == "26to49", nm.age.num] = 6
  dat.sub[dat.sub[, nm.age] == "50to69", nm.age.num] = 7
  dat.sub[dat.sub[, nm.age] == "70to79", nm.age.num] = 8
  dat.sub[dat.sub[, nm.age] == "80plus", nm.age.num] = 9
  dat.sub[, nm.age.num] = as.numeric(dat.sub[, nm.age.num])
  dat.sub[, "race_bin"] = ifelse(dat.sub[, nm.race] == "White", 
                                 1, 0)
  dat.sub
}