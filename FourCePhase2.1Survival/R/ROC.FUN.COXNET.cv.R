ROC.FUN.COXNET.cv=function(t0.all, data.t, data.v1, data.v2, data.v3, K=5){
  FPR0=seq(.01,.99,by=.01)
  yy0 = 0.5
  nn.t <- nrow(data.t)
  nn.v1 <- nrow(data.v1)
  nn.v2 <- nrow(data.v2)
  nn.v3 <- nrow(data.v3)
  #pp = ncol(data.t)
  pnum.t=data.t[,1]
  pnum.v1=data.v1[,1]
  pnum.v2=data.v2[,1]
  pnum.v3=data.v3[,1]
  set.seed(1234)
  pnum.v1.shuffle=sample(pnum.v1, replace=F)
  pnum.v2.shuffle=sample(pnum.v2, replace=F)
  pnum.v3.shuffle=sample(pnum.v3, replace=F)
  
  betahat = tryCatch(Est.ALASSO.GLMNET(data.t[,-1],fam0="Cox", w.b=NULL, Wi=NULL)$bhat.BIC,error=function(e) NA)
  if(is.na(betahat[1])!=1){
    ### v1
    nk.v1=floor(nn.v1/K)
    yyi.v1=array(NA, c(nn.v1, length(t0.all)))
    for(k in 1:K){
      pnum.v1.k = pnum.v1.shuffle[1:nk.v1 + (k-1)*nk.v1]
      pnum.t.k = setdiff(pnum.t,pnum.v1.k)
      beta.t = tryCatch(Est.ALASSO.GLMNET(data.t[data.t[,1]%in%pnum.t.k, -1],fam0="Cox", Wi=NULL)$bhat.BIC,error=function(e) NA)
      
      junk=do.call(cbind,lapply(t0.all, function(t0) 
        1-FUN.predict.cox(newZ=data.v1[data.v1[,1]%in%pnum.v1.k,-c(1:3)],beta=beta.t,X=data.t[,2],d=data.t[,3],Z=data.t[,-c(1:3)],t0=t0)
      ))
      yyi.v1[pnum.v1%in%pnum.v1.k,]=junk
    }
    junk.v1=lapply(t0.all, function(tt) ROC.Est.FUN(1*(data.v1[,3]*(data.v1[,2]<=tt)),yyi.v1[,tt],yy0,FPR0))
    junk.v1=do.call(cbind, junk.v1)
    yyi.v1=data.frame(patient_num=data.v1[,1], yyi.v1)
    
    ### v2
    nk.v2=floor(nn.v2/K)
    yyi.v2=array(NA, c(nn.v2, length(t0.all)))
    for(k in 1:K){
      pnum.v2.k = pnum.v2.shuffle[1:nk.v2 + (k-1)*nk.v2]
      pnum.t.k = setdiff(pnum.t,pnum.v2.k)
      beta.t = tryCatch(Est.ALASSO.GLMNET(data.t[data.t[,1]%in%pnum.t.k, -1],fam0="Cox", Wi=NULL)$bhat.BIC,error=function(e) NA)
      
      junk=do.call(cbind,lapply(t0.all, function(t0) 
        1-FUN.predict.cox(newZ=data.v2[data.v2[,1]%in%pnum.v2.k,-c(1:3)],beta=beta.t,X=data.t[,2],d=data.t[,3],Z=data.t[,-c(1:3)],t0=t0)
      ))
      yyi.v2[pnum.v2%in%pnum.v2.k,]=junk
    }
    junk.v2=lapply(t0.all, function(tt) ROC.Est.FUN(1*(data.v2[,3]*(data.v2[,2]<=tt)),yyi.v2[,tt],yy0,FPR0))
    junk.v2=do.call(cbind, junk.v2)
    yyi.v2=data.frame(patient_num=data.v2[,1], yyi.v2)
    
    ### v3
    nk.v3=floor(nn.v3/K)
    yyi.v3=array(NA, c(nn.v3, length(t0.all)))
    for(k in 1:K){
      pnum.v3.k = pnum.v3.shuffle[1:nk.v3 + (k-1)*nk.v3]
      pnum.t.k = setdiff(pnum.t,pnum.v3.k)
      beta.t = tryCatch(Est.ALASSO.GLMNET(data.t[data.t[,1]%in%pnum.t.k, -1],fam0="Cox", Wi=NULL)$bhat.BIC,error=function(e) NA)
      
      junk=do.call(cbind,lapply(t0.all, function(t0) 
        1-FUN.predict.cox(newZ=data.v3[data.v3[,1]%in%pnum.v3.k,-c(1:3)],beta=beta.t,X=data.t[,2],d=data.t[,3],Z=data.t[,-c(1:3)],t0=t0)
      ))
      yyi.v3[pnum.v3%in%pnum.v3.k,]=junk
    }
    junk.v3=lapply(t0.all, function(tt) ROC.Est.FUN(1*(data.v3[,3]*(data.v3[,2]<=tt)),yyi.v3[,tt],yy0,FPR0))
    junk.v3=do.call(cbind, junk.v3)
    yyi.v3=data.frame(patient_num=data.v3[,1], yyi.v3)
  }else{junk.v1=junk.v2=junk.v3=yyi.v1=yyi.v2=yyi.v3=NA}
  return(list(betahat=betahat, roc1=junk.v1, roc2=junk.v2, roc3=junk.v3, yyi.v1=yyi.v1, yyi.v2=yyi.v2, yyi.v3=yyi.v3))
}