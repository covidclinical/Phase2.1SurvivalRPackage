Est.Survfit.COXNET=function(dat.train, nm.event, multi.formulas, t0.all){
  
  junk=Est.ALASSO.GLMNET(dat.train[,-1],fam0="Cox", Wi=NULL)
  betahat=betahat.AIC = junk$bhat.AIC
  betahat.modBIC=junk$bhat.modBIC
  betahat.BIC=junk$bhat.BIC
  
  lamhat=junk$lambda.AIC
  score.app=do.call(cbind,lapply(t0.all, function(t0) 
    1-FUN.predict.cox(newZ=dat.train[,-c(1:3)],beta=betahat,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
  ))
  score.app=data.frame(patient_num=dat.train$patient_num, score.app)
  colnames(score.app)=c("patient_num", t0.all)
  return(list(betahat=betahat, betahat.modBIC=betahat.modBIC, betahat.AIC=betahat.AIC,betahat.BIC=betahat.BIC,lamhat=lamhat, score.app=score.app))
}