CV.Survfit.COXNET.sep=function(dat.train, dat.valid, betahat, lamhat, t0.all, nm.event, K, is.bt, yes.cv=T){
  pnum.train=dat.train[,1]
  pnum.valid=dat.valid[,1]
  nn=length(pnum.valid)
  nm.cov=names(betahat)[grepl("ns.calendar",names(betahat))!=1]

  X.train=dat.train[,"days_since_admission"]
  D.train=dat.train[,nm.event]
  Z.train= dat.train[,setdiff(colnames(dat.train), c("patient_num", "days_since_admission", nm.event))]
  
  X.valid=dat.valid[,"days_since_admission"]
  D.valid=dat.valid[,nm.event]
  Z.valid= dat.valid[,setdiff(colnames(dat.valid), c("patient_num", "days_since_admission", nm.event))]
  
  if(yes.cv==T){
    #lambda.grid=exp(seq(log(lamhat)-1, log(lamhat)+10,0.01))
    lambda.grid= 10^seq(-4,0,0.01)
    
  nk=floor(nn/K)
  yyi=array(NA, c(nn, length(t0.all)))
  yyi.sep.all=yyi.sep.cov=rep(NA,nn)
  
  for(k in 1:K){
    pnum.v = pnum.valid[1:nk + (k-1)*nk]
    if(k==K){pnum.v=pnum.valid[(nk+(k-2)*nk+1):length(pnum.valid)]}
    pnum.t = setdiff(pnum.train,pnum.v)
    junk=tryCatch(Est.ALASSO.GLMNET(dat.train[dat.train[,1]%in%pnum.t,-1], fam0="Cox", Wi=NULL, lambda.grid=lambda.grid),error=function(e) NA)
    beta.t = tryCatch(junk$bhat.modBIC,error=function(e) NA)
    junk=do.call(cbind,lapply(t0.all, function(t0) 
      1-FUN.predict.cox(newZ=dat.valid[dat.valid[,1]%in%pnum.v,-c(1:3)],beta=beta.t,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
    ))
    yyi[pnum.valid%in%pnum.v,]=junk
    yyi.sep.all[pnum.valid%in%pnum.v]=data.matrix(dat.valid[dat.valid[,1]%in%pnum.v,-c(1:3)])%*%beta.t
    yyi.sep.cov[pnum.valid%in%pnum.v]=data.matrix(dat.valid[dat.valid[,1]%in%pnum.v,nm.cov])%*%beta.t[nm.cov]
     }
  }else{
    yyi=do.call(cbind,lapply(t0.all, function(t0) 
      1-FUN.predict.cox(newZ=dat.valid[,-c(1:3)],beta=betahat,X=dat.train[,2],d=dat.train[,3],Z=dat.train[,-c(1:3)],t0=t0)
    ))
    yyi.sep.all=data.matrix(dat.valid[,-c(1:3)])%*%betahat
    yyi.sep.cov=data.matrix(dat.valid[,nm.cov])%*%betahat[nm.cov]
     }
  score.cv=data.frame(patient_num=pnum.valid, yyi)
  score.sep.all.cv=data.frame(patient_num=pnum.valid, yyi.sep.all)
  score.sep.cov.cv=data.frame(patient_num=pnum.valid, yyi.sep.cov)
 
  colnames(score.cv)=c("patient_num", t0.all)
  colnames(score.sep.all.cv)=
    colnames(score.sep.cov.cv)=
    c("patient_num", "score")
  if(is.bt==T){
    junk.bt=tryCatch(lapply(1:100, function(myseed){
    set.seed(myseed)
    dat.sample=dat.train[sample(1:dim(dat.train)[1], replace=T),]
    betahat.bt=tryCatch(Est.ALASSO.GLMNET(dat.sample[,-1], fam0="Cox", Wi=NULL, lambda.grid=lambda.grid)$bhat.modBIC,error=function (e) NA)
    if(length(betahat.bt)!=1){
    yyi.sep.cov=data.matrix(dat.sample[,nm.cov])%*%betahat.bt[nm.cov]
     
    score.sep.cov.bt=data.frame(patient_num=dat.sample$patient_num, score=yyi.sep.cov)
      }else{
      score.sep.cov.bt=NA}
    junk.bt=list(score.sep.cov.bt=score.sep.cov.bt, betahat.bt=betahat.bt)
    junk.bt
    }
    ), error=function(e) NA)
    }else{
      junk.bt=list(score.sep.cov.bt=NA, betahat.bt=NA)
      junk.bt
    }
  return(list(score.cv=score.cv, 
              score.sep.all.cv=score.sep.all.cv, 
              score.sep.cov.cv=score.sep.cov.cv,
              junk.bt=junk.bt
              #junk.bt=junk.bt#,
              #betahat.bt=betahat.bt,
              #beta.t=beta.t,
              #check=colnames(dat.sample)
              ))
}
