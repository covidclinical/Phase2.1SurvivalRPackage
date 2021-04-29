CV.Survfit.CSTAT.sep=function(dat.train, dat.valid, nm.event, K, is.bt, yes.cv=T){
  pnum.train=dat.train[,1]
  pnum.valid=dat.valid[,1]
  nn=length(pnum.valid)
  nm.cov=setdiff(colnames(dat.train),c("patient_num", "days_since_admission", nm.event))

  X.train=dat.train[,"days_since_admission"]
  D.train=dat.train[,nm.event]
  Z.train= dat.train[,setdiff(colnames(dat.train), c("patient_num", "days_since_admission", nm.event))]
  
  X.valid=dat.valid[,"days_since_admission"]
  D.valid=dat.valid[,nm.event]
  Z.valid= dat.valid[,setdiff(colnames(dat.valid), c("patient_num", "days_since_admission", nm.event))]
  
  if(yes.cv==T){
  cstat=NULL
  nn=dim(dat.valid)[1]
  nk=floor(nn/K)
  for(k in 1:K){
    pnum.v = pnum.valid[1:nk + (k-1)*nk]
    if(k==K){pnum.v=pnum.valid[(nk+(k-2)*nk+1):length(pnum.valid)]}
    pnum.t = setdiff(pnum.train,pnum.v)
    dat.train.k=dat.train[dat.train[,1]%in%pnum.t,-1]
    dat.valid.k=dat.valid[dat.valid[,1]%in%pnum.v,-1]
    
    junk=as.formula(paste0(
      paste0('Surv(days_since_admission,', nm.event,')~'),
      paste0(paste(colnames(Z.valid), collapse="+"))))
    
    train.fit <- coxph(junk,
                       x=TRUE, y=TRUE, method="breslow", data=dat.train.k)
    lpnew <- predict(train.fit, newdata=dat.valid.k) 
    Surv.rsp <- Surv(dat.train.k$days_since_admission, dat.train.k[,nm.event]) 
    Surv.rsp.new <- Surv(dat.valid.k$days_since_admission, dat.valid.k[,nm.event])
    cstat=c(cstat,UnoC(Surv.rsp, Surv.rsp.new, lpnew))
  }
  cstat.cv=mean(cstat,na.rm=T)
  }else{cstat.cv=NULL}
  
  if(is.bt==T){
    cstat.bt=unlist(lapply(1:100, function(myseed){
    set.seed(myseed)
    dat.bt=dat.train[sample(1:dim(dat.train)[1], replace=T),]
    train.bt <- coxph(junk,
                       x=TRUE, y=TRUE, method="breslow", data=dat.bt)
    lpnew.bt <- predict(train.bt, newdata=dat.bt) 
    Surv.rsp.bt <- Surv(dat.bt$days_since_admission, dat.bt[,nm.event]) 
    Surv.rsp.new.bt <- Surv(dat.bt$days_since_admission, dat.bt[,nm.event])
    UnoC(Surv.rsp.bt, Surv.rsp.new.bt, lpnew.bt)
    }))
  }else{cstat.bt=NA}
    
  return(list(cstat.cv=cstat.cv, cstat.bt=cstat.bt))
}
