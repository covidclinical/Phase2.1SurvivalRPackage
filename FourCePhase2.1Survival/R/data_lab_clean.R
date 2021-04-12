data_lab_clean=function (dat.x.raw, code.dict, nm.patient_num, nm.days_since_admission, 
                         nm.value, day) 
{
  dat.sub = dat.x.raw[dat.x.raw[, nm.days_since_admission] == 
                        day & dat.x.raw[, "concept_type"] %in% "LAB-LOINC", ]
  dat.sub$concept = paste(dat.sub$concept_type, dat.sub$concept_code, 
                          sep = ":")
  dat.sub = dat.sub[, c(nm.patient_num, "concept", nm.value)]
  #dat.sub.wide <- spread(dat.sub, concept, value)
  dat.sub.wide <- dat.sub %>%
    pivot_wider(values_from = value, names_from = concept, values_fn = mean)
  
  dat.sub.wide[dat.sub.wide == -99] = NA
  dat.sub.wide[dat.sub.wide == -999] = NA
  if ("LAB-LOINC:48065-7" %in% colnames(dat.sub.wide)) {
    dat.sub.wide[, "LAB-LOINC:48065-7:48066-5"] = 0.5 * dat.sub.wide[, 
                                                                     "LAB-LOINC:48065-7"]
  }
  if ("LAB-LOINC:48066-5" %in% colnames(dat.sub.wide)) {
    dat.sub.wide[, "LAB-LOINC:48065-7:48066-5"] = dat.sub.wide[, 
                                                               "LAB-LOINC:48066-5"]
  }
  dat.sub.wide = dat.sub.wide[, setdiff(colnames(dat.sub.wide), 
                                        c("LAB-LOINC:48065-7", "LAB-LOINC:48066-5"))]
  nm.lab = colnames(dat.sub.wide)[-1]
  code.dict.new = code.dict
  colnames(code.dict.new)[1] = "concept"
  code.dict.new[, 1] = paste("LAB-LOINC", code.dict.new[, 1], 
                             sep = ":")
  nm.lab = left_join(data.frame(concept = nm.lab), code.dict.new, 
                     by = "concept")
  nm.lab = nm.lab[, 2]
  nm.lab = gsub(" ", "_", nm.lab)
  nm.lab = gsub("\\(", "", nm.lab)
  nm.lab = gsub("\\)", "", nm.lab)
  nm.lab = gsub("-", "_", nm.lab)
  colnames(dat.sub.wide) = c(colnames(dat.sub.wide)[1], nm.lab)
  rownames(dat.sub.wide) = NULL
  data.frame(dat.sub.wide[, which(is.na(colnames(dat.sub.wide)) != 1)])
}