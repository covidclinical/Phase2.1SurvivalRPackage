fit.coxnet=function(t0.all, nm.event, multi.formulas,
                        dat.fit.t, dat.fit.v1, dat.fit.v2, dat.fit.v3){
  Xi.t=dat.fit.t$days_since_admission
  Di.t=dat.fit.t[,nm.event]
  
  Xi.v1=dat.fit.v1$days_since_admission
  Di.v1=dat.fit.v1[,nm.event]
  
  Xi.v2=dat.fit.v2$days_since_admission
  Di.v2=dat.fit.v2[,nm.event]
  
  Xi.v3=dat.fit.v3$days_since_admission
  Di.v3=dat.fit.v3[,nm.event]
  
    
    dat.tmp.t=dat.fit.t
    dat.tmp.v1=dat.fit.v1
    dat.tmp.v2=dat.fit.v2
    dat.tmp.v3=dat.fit.v3

  
    X.t<- model.matrix(multi.formulas, dat.tmp.t)[, -1]
    X.v1<- model.matrix(multi.formulas, dat.tmp.v1)[, -1]
    X.v2<- model.matrix(multi.formulas, dat.tmp.v2)[, -1]
    X.v3<- model.matrix(multi.formulas, dat.tmp.v3)[, -1]
    
    data.t=cbind(dat.tmp.t[,c(1:3)], X.t)
    data.v1=cbind(dat.tmp.v1[,c(1:3)], X.v1)
    data.v2=cbind(dat.tmp.v2[,c(1:3)], X.v2)
    data.v3=cbind(dat.tmp.v3[,c(1:3)], X.v3)
    
    junk=ROC.FUN.COXNET.cv(t0.all, data.t, data.v1, data.v2, data.v3)
    junk
}