lab_outlier_fun=function(dat.lab){
  nm.10lab=c("alanine_aminotransferase_ALT",     
             "albumin",
             "aspartate_aminotransferase_AST",
             "creatinine",
             "C_reactive_protein_CRP_Normal_Sensitivity",
             "total_bilirubin",
             "white_blood_cell_count_Leukocytes",
             "lymphocyte_count",
             "neutrophil_count", 
             "D_dimer")
  
  out.fun=function(x){
  mm=median(x,na.rm=T)
  UB=mm+10*mad(x, na.rm=T)
  ind.out=which(x>UB)
  x.new=x
  x.new[ind.out]=NA
  list(ind.out=ind.out, x.new=x.new)
  }
  
  nm.trans=c("alanine_aminotransferase_ALT","aspartate_aminotransferase_AST", 
             "C_reactive_protein_CRP_Normal_Sensitivity", "Ferritin", 
             "lactate_dehydrogenase_LDH","D_dimer")
  
  dat.lab.trans=dat.lab[,nm.10lab]
  dat.lab.trans=log(dat.lab.trans+0.5)

  dat.lab.trans.new=do.call(cbind,lapply(nm.10lab, function(ll) out.fun(dat.lab.trans[,ll])$x.new))
  colnames(dat.lab.trans.new)=nm.10lab
  dat.lab.new=exp(dat.lab.trans.new)-0.5
  ind.out=unique(unlist(lapply(nm.10lab, function(ll) out.fun(dat.lab.trans[,ll])$ind.out)))
  patient.num.outlier=dat.lab$patient_num[ind.out]
  dat.lab.new=data.frame(patient_num=dat.lab$patient_num, dat.lab.new)
  list(dat.lab.new=dat.lab.new, patient.num.outlier=patient.num.outlier)
}
