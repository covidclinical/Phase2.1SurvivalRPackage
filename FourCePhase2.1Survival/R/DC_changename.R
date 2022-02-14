change_nm <- function(dat){
  colnames(dat)=gsub("C_reactive_protein_CRP_Normal_Sensitivity", "CRP", colnames(dat))
  colnames(dat)=gsub("lymphocyte_count", "LYM", colnames(dat))
  colnames(dat)=gsub("D_dimer", "DD", colnames(dat))
  colnames(dat)=gsub("alanine_aminotransferase_ALT", "ALT", colnames(dat))
  colnames(dat)=gsub("aspartate_aminotransferase_AST", "AST", colnames(dat))
  colnames(dat)=gsub("total_bilirubin", "TB", colnames(dat))
  colnames(dat)=gsub("lactate_dehydrogenase_LDH", "LDH", colnames(dat))
  colnames(dat)=gsub("prothrombin_time_PT", "PT", colnames(dat))
  colnames(dat)=gsub("white_blood_cell_count_Leukocytes", "WBC", colnames(dat))
  colnames(dat)=gsub("age_group_new", "age", colnames(dat))
  colnames(dat)=gsub("race_new", "race", colnames(dat))
  colnames(dat)=gsub("cardiac_troponin_High_Sensitivity", "cardiac_high", colnames(dat))
  colnames(dat)=gsub("cardiac_troponin_Normal_Sensitivity", "cardiac_morm", colnames(dat))
  colnames(dat)=gsub("_CLASS", "", colnames(dat))
  colnames(dat)=gsub(".CLASS", "", colnames(dat))
  
  if("ALT"%in%colnames(dat) & "AST"%in%colnames(dat)){
    X.ALT=dat$ALT
    X.AST=dat$AST
    dat=data.frame(dat, AA=X.AST/X.ALT)
  }
  dat=dat[,grepl("ALT", colnames(dat))!=1]
  
  nm.lab.trans=c("CRP", "AST", "ALT", "DD")
  nm.med.trans=colnames(dat)[grepl("MED",colnames(dat))]
  dat[,colnames(dat)%in%c(nm.lab.trans, nm.med.trans)]=log(dat[,colnames(dat)%in%c(nm.lab.trans, nm.med.trans)]+0.5)
  dat0=dat
  dat0=dat0[dat0$days_since_admission!=0,]
  return(dat0)
}