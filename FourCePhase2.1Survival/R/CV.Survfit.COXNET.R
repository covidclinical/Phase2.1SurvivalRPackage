CV.Survfit.COXNET=function(dat.train, dat.valid, betahat, lamhat, t0.all, nm.event, K, yes.cv=T){
  pnum.train=dat.train[,1]
  pnum.valid=dat.valid[,1]
  nn=length(pnum.valid)
  
  X.train=dat.train[,"days_since_admission"]
  D.train=dat.train[,nm.event]
  Z.train= dat.train[,setdiff(colnames(dat.train), c("patient_num", "days_since_admission", nm.event))]
  
  X.valid=dat.valid[,"days_since_admission"]
  D.valid=dat.valid[,nm.event]
  Z.valid= dat.valid[,setdiff(colnames(dat.valid), c("patient_num", "days_since_admission", nm.event))]
  
  if(yes.cv==T){
    lambda.grid=exp(seq(log(lamhat)-1, log(lamhat)+1,0.01))
  nk=floor(nn/K)
  yyi=array(NA, c(nn, length(t0.all)))
  
  for(k in 1:K){
    pnum.v = pnum.valid[1:nk + (k-1)*nk]
    pnum.t = setdiff(pnum.train,pnum.v)
    beta.t = tryCatch(Est.ALASSO.GLMNET(dat.train[dat.train[,1]%in%pnum.t,-1], fam0="Cox", Wi=NULL, lambda.grid=lambda.grid)$bhat.AIC,error=function(e) NA)
    junk=do.call(cbind,lapply(t0.all, function(t0) 
      1-FUN.predict.cox(newZ=dat.valid[dat.valid[,1]%in%pnum.v,-c(1:3)],beta=beta.t,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
    ))
    yyi[pnum.valid%in%pnum.v,]=junk
  }
  }else{
    yyi=do.call(cbind,lapply(t0.all, function(t0) 
      1-FUN.predict.cox(newZ=dat.valid[,-c(1:3)],beta=beta.t,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
    ))}
  
  score.cv=data.frame(patient_num=pnum.valid, yyi)
  colnames(score.cv)=c("patient_num", t0.all)
  return(list(score.cv=score.cv))
}
