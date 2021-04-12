CV.Survfit.GLMNET=function(dat.train, dat.valid, betahat, lamhat, t0.all, nm.event, ipw=T, K, yes.cv=T){
  pnum.train=dat.train[,1]
  pnum.valid=dat.valid[,1]
  nn=length(pnum.valid)
  
  X.train=dat.train[,"days_since_admission"]
  D.train=dat.train[,nm.event]
  Z.train= dat.train[,setdiff(colnames(dat.train), c("patient_num", "days_since_admission", nm.event))]
  
  X.valid=dat.valid[,"days_since_admission"]
  D.valid=dat.valid[,nm.event]
  Z.valid= dat.valid[,setdiff(colnames(dat.valid), c("patient_num", "days_since_admission", nm.event))]
  
  score.cv=NULL
  for(tt in t0.all){
    Y.train=  I(X.train <= tt)*D.train
    dat.train.tt=cbind(pnum.train,Y.train, Z.train)
    
    Y.valid=  I(X.valid <= tt)*D.valid
    dat.valid.tt=cbind(pnum.valid, Y.valid, Z.valid)
    
  if(ipw==T){
    Gt.train=WGT.CEN(X.train, D.train, tt)
    Gt.valid=WGT.CEN(X.valid, D.valid, tt)
  }else{Gt.train=Gt.valid=NULL}
    
  if(yes.cv==T){
  lambda.grid=exp(seq(log(lamhat[names(lamhat)==tt])-1, log(lamhat[names(lamhat)==tt])+1,0.01))
  nk=floor(nn/K)
  yyi=rep(NA, nn)
    for(k in 1:K){
      pnum.v = pnum.valid[1:nk + (k-1)*nk]
      pnum.t = setdiff(pnum.train,pnum.v)
      beta.t = tryCatch(Est.ALASSO.GLMNET(dat.train.tt[dat.train.tt[,1]%in%pnum.t,-1],Wi=Gt.train[dat.train[,1]%in%pnum.t], lambda.grid=lambda.grid)$bhat.BIC,error=function(e) NA)
      yyi[pnum.valid%in%pnum.v] = g.logit(cbind(1,as.matrix(dat.valid.tt[dat.valid.tt[,1]%in%pnum.v,-(1:2)]))%*%beta.t)
    }
  }else{
  yyi=g.logit(cbind(1,as.matrix(dat.valid.tt[,-(1:2)]))%*%betahat[,colnames(betahat)==tt])
  }
  score.cv=cbind(score.cv, yyi)
  }
  score.cv=data.frame(patient_num=pnum.valid, score.cv)
  colnames(score.cv)=c("patient_num", t0.all)
  return(list(score.cv=score.cv))
}
