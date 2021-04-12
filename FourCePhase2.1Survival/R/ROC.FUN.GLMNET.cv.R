
ROC.FUN.GLMNET.cv  <- function(data.t, data.v1, data.v2, data.v3, 
                                wgti.t=NULL, wgti.v1=NULL, wgti.v2=NULL, wgti.v3=NULL, K=10){
  FPR0=seq(.01,.99,by=.01)
  yy0 = 0.5
  nn.t <- nrow(data.t); if(is.null(wgti.t)){wgti.t=rep(1,nn.t)}
  nn.v1 <- nrow(data.v1); if(is.null(wgti.v1)){wgti.v1=rep(1,nn.v1)}
  nn.v2 <- nrow(data.v2); if(is.null(wgti.v2)){wgti.v2=rep(1,nn.v2)}
  nn.v3 <- nrow(data.v3); if(is.null(wgti.v3)){wgti.v3=rep(1,nn.v3)}
  pp = ncol(data.t)
  pnum.t=data.t[,1]
  pnum.v1=data.v1[,1]
  pnum.v2=data.v2[,1]
  pnum.v3=data.v3[,1]
  set.seed(1234)
  pnum.v1.shuffle=sample(pnum.v1, replace=F)
  pnum.v2.shuffle=sample(pnum.v2, replace=F)
  pnum.v3.shuffle=sample(pnum.v3, replace=F)
  junk=tryCatch(Est.ALASSO.GLMNET(data.t[,-1],Wi=wgti.t))
  betahat = tryCatch(junk$bhat.BIC,error=function(e) NA)
  lamhat=tryCatch(junk$lambda.BIC, error=function(e) NA)
  if(is.na(betahat[1])!=1){
      lambda.grid.new=exp(seq(log(lamhat)-1, log(lamhat)+1,0.01))
      junk.v1=junk.v2=junk.v3=NULL
      ### v1
      nk.v1=floor(nn.v1/K)
      yyi.v1=rep(NA, nn.v1)
      for(k in 1:K){
      pnum.v1.k = pnum.v1.shuffle[1:nk.v1 + (k-1)*nk.v1]
      pnum.t.k = setdiff(pnum.t,pnum.v1.k)
      beta.t = tryCatch(Est.ALASSO.GLMNET(data.t[data.t[,1]%in%pnum.t.k,-1],Wi=wgti.t[data.t[,1]%in%pnum.t.k], lambda.grid=lambda.grid.new)$bhat.BIC,error=function(e) NA)
      yyi.v1[pnum.v1%in%pnum.v1.k] = g.logit(cbind(1,as.matrix(data.v1[pnum.v1%in%pnum.v1.k,-(1:2)]))%*%beta.t)
      }
      junk.v1=ROC.Est.FUN(data.v1[,2],yyi.v1,yy0,FPR0,wgti.v1)
      yyi.v1=data.frame(patient_num=data.v1[,1], score=yyi.v1)
      ### v2
      nk.v2=floor(nn.v2/K)
      yyi.v2=rep(NA, nn.v2)
      for(k in 1:K){
        pnum.v2.k = pnum.v2.shuffle[1:nk.v2 + (k-1)*nk.v2]
        pnum.t.k = setdiff(pnum.t,pnum.v2.k)
        beta.t = tryCatch(Est.ALASSO.GLMNET(data.t[data.t[,1]%in%pnum.t.k,-1],Wi=wgti.t[data.t[,1]%in%pnum.t.k], lambda.grid=lambda.grid.new)$bhat.BIC,error=function(e) NA)
        yyi.v2[pnum.v2%in%pnum.v2.k] = g.logit(cbind(1,as.matrix(data.v2[pnum.v2%in%pnum.v2.k,-(1:2)]))%*%beta.t)
      }
      junk.v2=ROC.Est.FUN(data.v2[,2],yyi.v2,yy0,FPR0,wgti.v2)
      yyi.v2=data.frame(patient_num=data.v2[,1], score=yyi.v2)
      
      ### for v3
      if(sum(pnum.t%in%pnum.v3)!=0){
      nk.v3=floor(nn.v3/K)
      yyi.v3=rep(NA, nn.v3)
      for(k in 1:K){
        pnum.v3.k = pnum.v3.shuffle[1:nk.v3 + (k-1)*nk.v3]
        pnum.t.k = setdiff(pnum.t,pnum.v3.k)
        beta.t = tryCatch(Est.ALASSO.GLMNET(data.t[data.t[,1]%in%pnum.t.k,-1],Wi=wgti.t[data.t[,1]%in%pnum.t.k], lambda.grid=lambda.grid.new)$bhat.BIC,error=function(e) NA)
        yyi.v3[pnum.v3%in%pnum.v3.k] = g.logit(cbind(1,as.matrix(data.v3[pnum.v3%in%pnum.v3.k,-(1:2)]))%*%beta.t)
      }
      junk.v3=ROC.Est.FUN(data.v3[,2],yyi.v3,yy0,FPR0,wgti.v3)
      }else{
      yyi.v3=g.logit(cbind(1,as.matrix(data.v3[,-(1:2)]))%*%betahat)
      junk.v3=ROC.Est.FUN(data.v3[,2],yyi.v3,yy0,FPR0,wgti.v3)
      }
      yyi.v3=data.frame(patient_num=data.v3[,1], score=yyi.v3)
  }else{junk.v1=junk.v2=junk.v3=yyi.v1=yyi.v2=yyi.v3=NA}
  return(list(betahat=betahat, roc1=junk.v1, roc2=junk.v2, roc3=junk.v3, 
              yyi.v1=yyi.v1, yyi.v2=yyi.v2, yyi.v3=yyi.v3))
}
