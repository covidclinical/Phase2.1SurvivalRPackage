data_lab_clean3=function (dat.x.raw, code.dict, nm.value, day) 
{
  dat.sub = dat.x.raw[dat.x.raw[, "days_since_admission"] == 
                        day & dat.x.raw[, "concept_type"] %in% "LAB-LOINC", ]
  dat.sub$concept = paste(dat.sub$concept_type, dat.sub$concept_code, 
                          sep = ":")
  dat.sub = dat.sub[, c("patient_num", "concept", nm.value)]
  dat.sub.wide <- dat.sub %>%
    pivot_wider(values_from = value, names_from = concept, values_fn = mean)
  #dat.sub.wide <- spread(dat.sub, concept, value)
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
  code.dict.new[,2]=as.character(code.dict.new[,2])
  code.dict.new=rbind(code.dict.new, c("LAB-LOINC:48065-7:48066-5", "D-dimer"))
  
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