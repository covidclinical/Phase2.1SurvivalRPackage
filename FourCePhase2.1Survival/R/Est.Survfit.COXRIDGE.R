Est.Survfit.COXRIDGE=function(dat.train, nm.event, multi.formulas, t0.all){
  
  junk=Est.ALASSO.GLMRIDGE(dat.train[,-1],fam0="Cox", Wi=NULL)
  betahat= junk$bhat

  score.app=do.call(cbind,lapply(t0.all, function(t0) 
    1-FUN.predict.cox(newZ=dat.train[,-c(1:3)],beta=betahat,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
  ))
  score.app=data.frame(patient_num=dat.train$patient_num, score.app)
  colnames(score.app)=c("patient_num", t0.all)
  return(list(betahat=betahat, score.app=score.app))
}