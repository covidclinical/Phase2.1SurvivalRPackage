Est.Survfit.coefficient.GLMNET=function(dat.train, betahat,nm.event,  t0.all){
  X.train=dat.train[,"days_since_admission"]
  D.train=dat.train[,nm.event]
  Z.train= dat.train[,setdiff(colnames(dat.train), c("patient_num", "days_since_admission", nm.event))]
  
  score.app=NULL
    for(tt in t0.all){
    print(tt)
    Y.train=  I(X.train <= tt)*D.train
    dat.train.tt=cbind(Y.train, Z.train)
    score.app=cbind(score.app,g.logit(cbind(as.matrix(dat.train.tt[,-1]))%*%betahat))
    }
  score.app=data.frame(patient_num=dat.train[,"patient_num"], score.app)
  colnames(score.app)=c("patient_num", t0.all)
  return(list(score.app=score.app))
}