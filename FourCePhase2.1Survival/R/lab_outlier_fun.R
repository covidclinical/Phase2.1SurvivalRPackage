lab_outlier_fun=function(dir.input){
  dat.x.raw = read.csv(paste0(dir.input, "/LocalPatientObservations.csv"))
  nm.patient_num = "patient_num"
  nm.days_since_admission = "days_since_admission"
  nm.value = "value"
  dat.check=data_lab_clean2(dat.x.raw, nm.patient_num, nm.days_since_admission, 
                            nm.value, day=c(0:1)) 
  
  out.fun=function(x){
    mm=median(x,na.rm=T)
    UB=mm+50* mad(x, na.rm=T)
    #ind.out=which((x>UB)|x>10000)
    ind.out=which((x>100*mm)|x>10000)
    #ind.out=which(x>10000)
    ind.out
  }
  
  nm.lab.all=colnames(dat.check)[-1]
  ind.out=unique(unlist(lapply(nm.lab.all, function(ll) out.fun(dat.check[,ll]))))
  patient_num.out=unique(dat.check$patient_num[ind.out])
  patient_num.out
}
