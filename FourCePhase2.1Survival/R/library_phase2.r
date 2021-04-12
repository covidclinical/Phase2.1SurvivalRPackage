g.logit = function(xx){exp(xx)/(exp(xx)+1)}
logit = function(xx){log(xx/(1-xx))}
dg.logit = function(xx){exp(xx)/(exp(xx)+1)^2}


logitlik.fun = function(bet.mat,dat){
  yi = dat[,1]; xi = dat[,-1]; pi.mat = g.logit(cbind(1,xi)%*%bet.mat) ## N x B
  apply(log(pi.mat)*yi + log(1-pi.mat)*(1-yi),2,sum)
}

WGT.CEN <- function(Ti, Di,t0)
{
  ## ================================================================ ##
  ## ============== KM Estimator of Censoring Survival ============== ##
  ## ================================================================ ##
  Ghat.FUN <- function(tt, Ti, Di,type='fl')
  {
    surv.fit <- survfit(Surv(Ti,1-Di)~1,se.fit=F,type=type)
    surv.ti <- surv.fit$time; surv.fit.surv = surv.fit$surv
    surv.til <- sort(surv.fit$time); surv.fit.surv = surv.fit.surv[order(surv.ti)]
    tmpind <- sum.I(tt,">=",surv.til) + 1
    c(1,surv.fit.surv)[tmpind]
  }
  
  N <- length(Ti)
  Ghat.Ti <- Ghat.FUN(Ti,Ti,Di)
  Ghat.t0 <- Ghat.FUN(t0,Ti,Di)
  Wi <- rep(0,length(Ti))
  Wi[Ti <= t0] <- Di[Ti<=t0]/Ghat.Ti[Ti<=t0]
  Wi[Ti >  t0] <- 1/Ghat.t0
  Wi
}

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}


## ================================================================================================================== ##
## estimate beta w/ adaptive LASSO regularization (if regularize=T) or standard logistic regression (if regularize=F) ##
## data: 1st column y; remaining x; nopen.ind indexes which subset of x should not be penalized                       ##
## Wi: weights for resampling if interested in obtaining standard errors                                              ## 
## ================================================================================================================== ##
Est.ALASSO.GLM = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize,pwr=0.1){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
    lam.ridge = pp/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
    bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family="binomial"))); ##print(bini)
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    ## glmpath provides solution path for a range of penalty parameters ##
    tmpfit = glmpath(x.t,y,nopenalty.subset=nopen.ind,family=binomial,weight=Wi,standardize=F,min.lambda=0,lambda2=lam.ridge)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda")
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(predict(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+min(sum(y)^pwr,log(sum(y)))*df.all 
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=binomial,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(out)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}


Est.ALASSO.GLM2 = function(data,Wi=NULL,rtn="EST",nopen.ind=NULL,regularize=yes.regularize, pwr=0.1, adaptive = TRUE, alter_initial = FALSE){
  data = as.matrix(data); y = data[,1]; x = data[,-1,drop=F]; nn=length(y); if(is.null(Wi)){Wi=rep(1,nn)}; pp = ncol(x)
  if(regularize){
    if(!alter_initial) {
      ## adaptive lasso (aLASSO) with initial estimator obtained via ridge (bini); w.b creates adaptive weights for aLASSO ##
      lam.ridge = log(pp)/nn; ##bini = plr(x,y,lambda=lam.ridge,weights=Wi)$coef; 
      bini = as.vector(coef(glmnet(x,y,weights=Wi,alpha=0,standardize=F,lambda=lam.ridge,family="binomial"))); ##print(bini)
      
    }
    else {
      bini = glmnet(x,y,weights=Wi,family="binomial",alpha=0,offset=offset)
      lam.xx=svd(t(x)%*%x/nn)$d; tmpdf = apply(lam.xx/(lam.xx+VTM(bini$lambda,pp)),2,sum)
      tmpind = which.min(deviance(bini)+2*tmpdf); 
      lam.ridge = bini$lambda[tmpind]
    }
    
    w.b = 1/abs(bini[-1]); x.t = x/VTM(w.b,nrow(x))
    
    if(!adaptive) {
      w.b = rep(1, ncol(x)); lam.ridge = 0;
    }
    
    ## glmpath provides solution path for a range of penalty parameters ##
    ## standarization 
    tmpfit = glmpath(x.t,y,nopenalty.subset=nopen.ind,family=binomial,weight=Wi,standardize=F,min.lambda=0,lambda2=lam.ridge)
    lam.all = c(seq(min(tmpfit$lambda),max(tmpfit$lambda),length=500))
    b.all = predict(tmpfit, s=lam.all, type="coefficients",mode="lambda")
    b.all = b.all/VTM(c(1,w.b),nrow(b.all)); m0 = length(lam.all)
    ## ================================================================================ ##
    ## calculates degree of freedom for all beta's (corresponding to different lam.all) ##
    ## ================================================================================ ##
    df.all = apply(b.all[,-1,drop=F]!=0,1,sum); x = as.matrix(x)
    ## =============================================================================================================== ##
    ## calculates modified BIC, log(n) is modified as min{sum(y)^0.1, log(n)} to avoid over shrinkage in finite sample ##
    ## =============================================================================================================== ##
    BIC.lam = -2*apply(predict(tmpfit,newx=x.t,newy=y,s=lam.all,type="loglik",mode="lambda"),2,sum)+min(sum(y)^pwr,log(sum(y)))*df.all 
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]; lamhat = lam.all[m.opt]
  }else{
    ## ========================================================================= ##
    ## if regularize = F, then use standard logistic regression w/o penalization ##
    ## ========================================================================= ##
    bhat=bini=glm(y~x,family=binomial,weight=Wi)$coef;lamhat = 0; lam.all=BIC.lam=b.all=NULL
  }
  out = c("b"=bhat, "bini"=bini,"lamhat"=lamhat,"lam.all"=lam.all,"BIC.lam"=BIC.lam,"b.all"=b.all)
  if(rtn=="EST"){return(out)}else{return(list(out,"b.all"=b.all,"lam.all"=lam.all,"fit"=tmpfit,"BIC.lam"=BIC.lam))}
}


SIM.FUN = function(nn,rtn="data.t")
  {
	xx = mvrnorm(nn,mu=rep(0,p.x),Sigma=Sig0.X); 
	icd.B = rbinom(nn,size=1,prob=g.logit(xx[,1]*3+2))
	xx[,1] = (xx[,1]*(xx[,1]>0)+(rexp(nn,rate=0.1)+5)*rbinom(nn,size=1,prob=0.1))*icd.B
	prob.x = g.logit(-alp0+c(xx%*%beta0)-3*(1-icd.B)+0.2*(xx[,1]>15))
	yy = rbinom(nn,prob=prob.x,size=1); dat = cbind(yy,xx)
	if(rtn=="data.t"){return(dat)}else{
		zz = rbinom(nn, size=2, prob=g.logit(log(maf)+gam.z*yy))
		return(list(dat,cbind("D"=yy,"P.x"=prob.x,"G"=zz)))}
  }


ROC.Est.FUN <- function(Di,yyi,yy0,fpr0=NULL,wgti=NULL,yes.smooth=F)
  {
    out.yy <- out.pp <- out.AUC <- out.TPR <- out.FPR <- out.PPV <- out.NPV <- NULL
    if(is.null(wgti)){wgti=rep(1,length(Di))}; yyi = as.matrix(yyi); pp=ncol(as.matrix(yyi));  
    mu0 = sum(wgti*(1-Di))/sum(wgti); mu1 = 1-mu0  
    for(k in 1:pp)
      {
       yy = yy0; 
       if(!is.null(fpr0)){
         tpr.all = S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
         fpr.all = S.FUN(yyi[,k],Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth);
         TPR = approx(c(0,fpr.all,1),c(0,tpr.all,1),fpr0,method="linear",rule=2)$y; 
         TPR = c(S.FUN(yy0,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth), TPR); 
          yy = c(yy,Sinv.FUN(fpr0,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth))           
         FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
       }else{
         TPR = S.FUN(yy,Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth); 
         FPR = S.FUN(yy,Yi=yyi[,k],(1-Di)*wgti,yes.smooth=yes.smooth)
       }
       out.yy = cbind(out.yy, yy)
       out.pp = cbind(out.pp, S.FUN(yy,Yi=yyi[,k],wgti,yes.smooth=yes.smooth))
       out.TPR = cbind(out.TPR,  TPR);  out.FPR  <- cbind(out.FPR,  FPR)
       PPV <- 1/(1+FPR*mu0/(TPR*mu1)); NPV <- 1/(1+(1-TPR)*mu1/((1-FPR)*mu0))
       out.PPV <- cbind(out.PPV, PPV); out.NPV <- cbind(out.NPV, NPV)
       #AUC <- sum((sum.I(yyi[,k],"<=",Yi=yyi[,k],Vi=Di*wgti)+sum.I(yyi[,k],"<",Yi=yyi[,k],Vi=Di*wgti))*(1-Di)*wgti/2
       #             )/(sum((1-Di)*wgti)*sum(Di*wgti))
       AUC = sum(S.FUN(yyi[,k],Yi=yyi[,k],Di*wgti,yes.smooth=yes.smooth)*(1-Di)*wgti)/sum((1-Di)*wgti)
       out.AUC <- c(out.AUC, AUC)
     }
    out = c(out.AUC,out.yy,out.pp,out.FPR,out.TPR,out.PPV,out.NPV)
    out
  }


S.FUN <- function(yy,Yi,Di,yes.smooth=F)
  {
  	if(yes.smooth){
		Y1i = Yi[Di==1]; n1 = sum(Di); bw = bw.nrd(Y1i)/n1^0.6
		c(t(rep(1/n1,n1))%*%pnorm((Y1i-VTM(yy,n1))/bw))
  	}else{
		return((sum.I(yy,"<",Yi,Vi=Di)+sum.I(yy,"<=",Yi,Vi=Di))/sum(Di)/2)
  	}
    ##sum.I(yy,"<=",Yi,Vi=Di)/sum(Di)
  }

Sinv.FUN <- function(uu,Yi,Di,yes.smooth=F)
  {
    yy0<-unique(sort(Yi,decreasing=T)); ss0 <- S.FUN(yy0,Yi,Di,yes.smooth=yes.smooth) 
    return(approx(ss0[!duplicated(ss0)],yy0[!duplicated(ss0)],uu,method="linear",f=0,rule=2)$y)
  }



VTM<-function(vc, dm){
    matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)
## sum_i I(yy FUN Yi)Vi
# Vi weight
  {
    if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
    # for each distinct ordered failure time t[j], number of Xi < t[j]
    pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')    
    if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos # number of Xi>= t[j]
    if (!is.null(Vi)) {
	   ## if FUN contains '=', tmpind is the order of decending
        if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
        ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
        Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
        return(rbind(0,Vi)[pos+1,])
    } else return(pos)
  }



