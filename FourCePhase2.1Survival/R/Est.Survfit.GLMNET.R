Est.Survfit.GLMNET=function(dat.train, nm.event, multi.formulas, t0.all, ipw=T){
  X.train=dat.train[,"days_since_admission"]
  D.train=dat.train[,nm.event]
  Z.train= dat.train[,setdiff(colnames(dat.train), c("patient_num", "days_since_admission", nm.event))]
  
  betahat=lamhat=score.app=NULL
    for(tt in t0.all){
    print(tt)
    Y.train=  I(X.train <= tt)*D.train
    dat.train.tt=cbind(Y.train, Z.train)
    if(ipw==T){Gt.train=WGT.CEN(X.train, D.train, tt)}else{Gt.t=NULL}
    junk=tryCatch(Est.ALASSO.GLMNET(dat.train.tt,Wi=Gt.train), error=function(e) NA)
    score.app=cbind(score.app,g.logit(cbind(1,as.matrix(dat.train.tt[,-1]))%*%junk$bhat.AIC))
    betahat = cbind(betahat,tryCatch(junk$bhat.AIC,error=function(e) NA))
    lamhat=c(lamhat, tryCatch(junk$lambda.AIC, error=function(e) NA))
    }
  score.app=data.frame(patient_num=dat.train[,"patient_num"], score.app)
  colnames(betahat)=t0.all
  names(lamhat)=t0.all
  colnames(score.app)=c("patient_num", t0.all)
  return(list(betahat=betahat, lamhat=lamhat, score.app=score.app))
}