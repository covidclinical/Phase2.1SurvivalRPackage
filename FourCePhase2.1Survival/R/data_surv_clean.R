data_surv_clean=function (dat.surv.raw, nm.patient_num, nm.days_since_admission, 
          daymax, nm.event) 
{
  dat.sub = dat.surv.raw[dat.surv.raw[, nm.days_since_admission] <= 
                           daymax, ]
  patient_num.list = sort(unique(dat.sub[, nm.patient_num]))
  data.surv.clean.kern = function(dat.sub, patient_num) {
    dat.tmp = dat.sub[dat.sub[, nm.patient_num] == patient_num, 
    ]
    dat.event.tmp = dat.tmp[which(dat.tmp[, nm.event] == 
                                    1), ]
    if (dim(dat.event.tmp)[1] != 0) {
      dat.end.tmp = dat.event.tmp[which.min(dat.event.tmp[, 
                                                          nm.days_since_admission]), c(nm.patient_num, 
                                                                                       nm.days_since_admission, nm.event)]
    }
    else {
      dat.end.tmp = dat.tmp[which.max(dat.tmp[, nm.days_since_admission]), 
                            c(nm.patient_num, nm.days_since_admission, nm.event)]
    }
    dat.end.tmp
  }
  dat.res = do.call(rbind, lapply(patient_num.list, function(xx) data.surv.clean.kern(dat.sub, 
                                                                                      xx)))
  rownames(dat.res) = NULL
  dat.res
}