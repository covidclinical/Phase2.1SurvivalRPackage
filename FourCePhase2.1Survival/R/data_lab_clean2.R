data_lab_clean2=function (code.dict, dat.x.raw, nm.value, day) 
{
  
  code.dict = apply(code.dict, 2, as.character)
  combine.set = c("48065-7", "48066-5")
  combine.nm = paste(combine.set, collapse = ":")
  code.dict = data.frame(rbind(code.dict, c(combine.nm, "D-dimer")))
  
  
  
  
  dat.sub = dat.x.raw[dat.x.raw[, "days_since_admission"] %in%
                        day & dat.x.raw[, "concept_type"] %in% "LAB-LOINC", ]
  dat.sub$concept = paste(dat.sub$concept_type, dat.sub$concept_code, 
                          sep = ":")
  dat.sub = dat.sub[, c("patient_num", "concept", nm.value)]
  dat.sub=data.frame(data.table(dat.sub)[, lapply(.SD, mean, na.rm=T), by=c("patient_num", "concept")])
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
  nm.lab = left_join(data.frame(concept = nm.lab), code.dict.new, 
                     by = "concept")
  nm.lab = nm.lab[, 2]
  nm.lab = gsub(" ", "_", nm.lab)
  nm.lab = gsub("\\(", "", nm.lab)
  nm.lab = gsub("\\)", "", nm.lab)
  nm.lab = gsub("-", "_", nm.lab)
  colnames(dat.sub.wide) = c(colnames(dat.sub.wide)[1], nm.lab)
  rownames(dat.sub.wide) = NULL
  dat.sub.wide=dat.sub.wide[, which(is.na(colnames(dat.sub.wide)) != 1)]
  data.frame(dat.sub.wide)
}