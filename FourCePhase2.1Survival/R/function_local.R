simul <- function(N, X, betas, lambda, rho, beta, rateC,dist_family='weibull'){
  ### get Weibull/Gompertz distribution with cumulative baseline hazard: H(t)=lambda*t^rho/H(t)=lambda/rho*(exp(rho*t)-1)
  ### rho changes across strata, for simplicity only assume 2 strata here
  ### formula to simulate time with Weibull, T using uniform distribution, v: (-log(v)/(lambda*exp(beta*X)))^(1/rho)
  ### to simulate time with Gompertz: ((1/rho)*log(1-(rho*log(v))/(lambda*exp(beta*X))))
  Tlat = rep(0,N)
  cbh = rep(0,N)
  v <- runif(n=N)
  # Weibull latent time
  rho0 = rho
  if(dist_family=='weibull'){
    Tlat = (-log(v)/(lambda*exp(as.matrix(X)%*%as.matrix(beta))))^(1/rho0)
  }
  if(dist_family=='gompertz'){
    Tlat = (1/rho0)*log(1-(rho0*log(v))/(lambda*exp(as.matrix(X)%*%as.matrix(beta))))
  }
  # censoring times
  #C1 <- rexp(n=N, rate=rateC)
  C2 = runif(n=N,quantile(Tlat,0.5),quantile(Tlat,0.9))
  #C2 = runif(n=N,3,6)
  # follow-up times and event indicators
  #time <- pmin(Tlat, C1, C2)
  #status <- as.numeric(Tlat <= C1 & Tlat <= C2)
  time <- pmin(Tlat,C2)
  status <- as.numeric(Tlat <= C2)
  # get theoretical value of cumulative baseline hazard and baseline hazard at each time
  return(cbind(data.frame(id=1:N,y=time,failed=status),X))
}


cbh = function(t,bhest){
  ### cumulative baseline hazard evaluated at a given time, based on output of basehaze function in R
  y = bhest$hazard[findInterval(t,bhest$time)]
  return(y)
}

ar1_cor <- function(n, rho) {
  ### construct autoregressive correlation structure
  exponent <- abs(matrix(1:n-1,nrow=n,ncol=n,byrow=TRUE)-(1:n-1))
  rho^exponent
}

fw = function(w,fw_COEF,fw_INTERCEPT){
  ### shape function
  fw_COEF*cos(w*pi/300)+fw_INTERCEPT
}
fw1 = function(w){
  ### interacting part of shape function
  sin(w*pi/300)
}

generate_data = function(n_patient,n_clinical_feature,betas,T0,
                         missing_rate,lambda,rho,rateC,Xcov,fw_ceof,fw_intercept,
                         d_phi=3,dist_family='weilbull'){
  ### the main function used to generate survival data
  ### Input are parameters 
  ### Output: dat.true: real but unobserved covariates used to generate survival time; dat.obs: observable covariates
  X0 = mvrnorm(n_patient, mu=rep(0,n_clinical_feature),Sigma=Xcov) # mvn with cov=ar(1)
  colnames(X0) = paste0('X',1:ncol(X0))
  miss_ind = matrix(0,nrow=nrow(X0),ncol=ncol(X0)) # missing indicator
  colnames(miss_ind) = paste0('M',1:ncol(X0))
  calendar_time = sample(1:T0, n_patient, replace = T) # sample calendar time
  myfw = data.frame(fw=fw(calendar_time,fw_coef,fw_intercept)) # generate shape function
  myfw1 = fw1(calendar_time)
  X = cbind(X0,miss_ind) # combine startum, mvn & missing indicator
  for (j in 1:n_clinical_feature){
    miss_indx = X0[,j]<1
    X[miss_indx,paste0('M',j)]=rbinom(sum(miss_indx),1,missing_rate[j])
    X[X[,n_clinical_feature+j]==1,paste0('X',j)]=0
  }
  interact = myfw1*X # add interaction term
  colnames(interact) = paste0('fw1.',colnames(X))
  phi = data.frame(ns(calendar_time,df=d_phi)) # add basis for calendar time
  colnames(phi) = paste0('phi',1:ncol(phi))
  # interaction with covariate
  X_ob = X
  for (j in 1:d_phi) {
    inter_basis = phi[,j]*X
    colnames(inter_basis) = paste0(colnames(phi)[j],".",colnames(X))
    X_ob = cbind(X_ob,inter_basis)
  }
  # interaction with missing indicator
  X_inter = cbind(X,interact)
  X_true = cbind(X_inter,myfw,myfw1) # combine coviate, missing indicator, interaction and shape fucntion
  mybetas = matrix(0,ncol=ncol(X_true))
  colnames(mybetas) = colnames(X_true)
  for (p in colnames(betas)){
    mybetas[1,p]=betas[1,p]
  }
  beta.X = X_inter%*%mybetas[,colnames(X_inter)]
  # use the weibull distribution to get survival time and status
  dat <- simul(N=n_patient, lambda=lambda, rho=rho, 
               beta=as.vector(mybetas), X=X_true, 
               rateC=rateC,dist_family=dist_family)
  dat.true = cbind(dat,calendar_time)
  dat.obs = cbind(subset(dat,select=c(id,y,failed)), X_ob, phi, calendar_time)
  cbh=get_true_cbh(dat.true$y,lambda,rho,dist_family=dist_family)$hazard
  return(list(dat.true=dat.true,dat.obs=dat.obs,cbh=cbh,beta.X=beta.X,fw=betas$fw*myfw))
}

getcbh = function(mymodel){
  ### get cumulative baseline hazard from adaptive lasso fit
  basesurvival = as.data.frame(basesurv(mymodel,center=FALSE))
  cbhest = data.frame(time=basesurvival$time)
  cbhest = cbind(cbhest,-log(subset(basesurvival,select = c(-time))))
  colnames(cbhest)=gsub('survival','bch',colnames(cbhest))
  return(cbhest)
}

log_lkh = function(eta,xi,beta,dat.obs,psi,d_phi,knot){
  ### compute full log likelihood of cox model
  psi = t(psi)
  cbhhat.obj = cbhhat(dat.obs$y,eta,knot)
  phi = t(dplyr::select(dat.obs,'phi1':paste0('phi',d_phi)))
  nonx_col=c('y','failed','calendar_time',paste0('phi',1:d_phi))
  X = t(dplyr::select(dat.obs,-nonx_col))
  lglkh = -mean(dat.obs$failed*(t(eta)%*%psi+t(xi)%*%phi+t(beta)%*%X)-
                  cbhhat.obj$int*exp(t(xi)%*%phi+t(beta)%*%X))
  return(lglkh)
}
grr = function(eta,xi,beta,dat.obs,psi,d_phi,knot){
  ### auxiliary function for optim, which compute the gradient of our -log-likelihood
  d_psi = length(eta)
  cbhhat.obj = cbhhat(dat.obs$y,eta,knot)
  phi = t(dplyr::select(dat.obs,'phi1':paste0('phi',d_phi)))
  nonx_col=c('y','failed','calendar_time',paste0('phi',1:d_phi))
  X = t(dplyr::select(dat.obs,-nonx_col))
  psi = t(psi)
  Delta = dat.obs$failed
  Z = as.numeric(t(xi)%*%phi+t(beta)%*%X)
  grad.lambda = matrix(0,d_psi,length(dat.obs$y))
  for (d in 1:d_psi){
    grad.lambda[d,] = (dat.obs$y/2)*(cbhhat.obj$weights%*%(cbhhat.obj$bh.eva*cbhhat.obj$psi.eva[,d,]))#colSums(cbhhat.obj$psi.eva[,d,]*(cbhhat.obj$bh.eva*cbhhat.obj$weights*dat.obs$y/2))
  }
  grad = -rowMeans(sweep(psi,MARGIN=2,Delta,'*')-sweep(grad.lambda,exp(Z),MARGIN=2,'*'))
  return(grad)
}

extr_cbh = function(mycbh,stratum){
  ### extract cumulative baseline hazard for a stratum
  output = data.frame(time=mycbh$time,"hazard"=mycbh[paste0('bch.S.',stratum)])
  colnames(output) = c('time','hazard')
  return(output)
}
get_true_cbh = function(time,lambda,rho,dist_family='weibull'){
  ### given time and degree for weibull dsitrbution
  ### output real cumulative baseline hazard
  if (dist_family == 'weibull'){
    hazard = lambda*(time)^(rho)
  }
  if (dist_family == 'gompertz'){
    hazard = lambda/rho*(exp(rho*time)-1)
  }
  return(data.frame(time=time,hazard=hazard))
}

range_scale = function(x){
  return((x-min(x))/(max(x)-min(x)))
}
psi_eval = function(t,time){
  return((t-min(time))/(max(time)-min(time)))
}
const_psi = function(cbhaz,time,n_node,scale=F){
  ### input: cumulative baseline hazard, time want to evaluate and dimension of Talor polynormial
  ### scale time and return Taylor polynormial of time  
  output = data.frame('time'=time,'hazard'=cbh(time,cbhaz))
  cutoff = seq(0,1,length.out=n_node+2)
  knot = quantile(time,cutoff)
  knot[length(knot)] = Inf
  # psi0 is the interecept term
  Psi = com_psi(time,knot)
  return(list(cbh=output,psi=Psi,knot=knot))
}

const_plot = function(dat1,dat2,dat3=NULL,scale=0,
                      label1='pred',label2='true',label3=NULL,xlab="",ylab="",title=""){
  ### input: pred/true: data.frame containing two columns, time and hazard
  ### plot predicted and true value
  if(length(dat3)==0){
    g=ggplot()+
      geom_point(data=dat1,aes(x=time,y=hazard,color=label1),alpha=0.5)+
      geom_point(data=dat2,aes(x=time,y=hazard*exp(scale),color=label2),alpha=0.5)+
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", values = c("red","blue"))+
      ggtitle(title)
  }else{
    g=ggplot()+
      geom_point(data=dat1,aes(x=time,y=hazard,color=label1),alpha=0.5)+
      geom_point(data=dat2,aes(x=time,y=hazard*exp(scale),color=label2),alpha=0.5)+
      geom_point(data=dat3,aes(x=time,y=hazard*exp(scale),color=label3),alpha=0.5)+
      xlab(xlab)+
      ylab(ylab)+
      scale_color_manual(name="", values = c("red","blue","green"))+
      ggtitle(title)
  }
  return(g)
}

get_local_est = function(dat.obs,num_bh_knot,d_phi){
  ### get estimation of parameters at local sites
  ### input: dat.obs: the covariate matrix; 
  ###        num_bh_knot: number of nodes used in piecewise-exponential function to approximate the baseline hazard function 
  ###        d_phi: dimensions for phi, the degrees of freedom of the cubic splines of calendar time 
  ### output: all estimations needed to compute summary statistics
  X_obs = dat.obs
  nonx_col=c('y','failed','calendar_time')
  xdata = dplyr::select(X_obs,-nonx_col)
  # fit stratified adaptive lasso
  penalty0 = data.frame(matrix(1,1,ncol(xdata)))
  colnames(penalty0)=colnames(xdata)
  for (i in 1:d_phi) {penalty0[paste0('phi',i)]=0}
  # first fit a ridge regression
  ridge_model = cv.glmnet(x=as.matrix(xdata),
                          y=Surv(X_obs$y, X_obs$failed),
                          family = 'cox',alpha=0,
                          penalty.factor = penalty0)
  # adaptivelasso, use absolute value coefficients of ridge as reciprocal of lasso penalty
  # alasso_model = penalized(unpen_part,
  #                          xdata,
  #                          data=X_obs,lambda1=as.vector(1/abs(coefficients(ridge_model)[-(1:d_phi)])),lambda2=0)
  penalty1 = as.vector(penalty0) / as.vector(abs(coef(ridge_model,s='lambda.min'))) 
  alasso_model = cv.glmnet(x=as.matrix(xdata),
                           y=Surv(X_obs$y, X_obs$failed),
                           family = 'cox',alpha=1,
                           penalty.factor = as.vector(penalty1))
  mymodel = as.formula(paste0('Surv(y, failed) ~',paste0(colnames(xdata),collapse='+')))
  cox_model <- coxph(mymodel, data=X_obs)
  cox_model$coefficients=coef(alasso_model,s='lambda.min')
  # extract beta
  betahat = as.matrix(coef(alasso_model,s='lambda.min')[setdiff(colnames(xdata),c(paste0('phi',1:d_phi))),])
  phi = as.matrix(dplyr::select(X_obs,paste0('phi',1:d_phi)))
  xi = as.matrix(coef(alasso_model,s='lambda.min')[paste0('phi',1:d_phi),])
  cumbh = basehaz(cox_model,centered = FALSE)
  psi.obj = const_psi(cumbh,X_obs$y,num_bh_knot)
  #mycbh = psi.obj$cbh
  psi = psi.obj$psi
  psi.knot = psi.obj$knot
  eta_initial = rep(.1,ncol(psi))
  eta = optim(eta_initial,fn=log_lkh,xi=xi,beta=betahat,dat.obs=X_obs,
              psi=psi,d_phi=d_phi,knot=psi.knot,gr=grr,method='L-BFGS-B',lower=c(rep(-Inf,length(eta_initial))))$par
  cbhhat.obj = cbhhat(X_obs$y,eta,psi.knot)
  myX = as.matrix(dplyr::select(xdata,-paste0('phi',1:d_phi)))
  return(list(X=t(myX),beta=betahat,eta=eta,xi=xi,phi=t(phi),
              psi=t(psi),
              alasso = alasso_model,cbh=cumbh,
              psi.obj=psi.obj,
              cbhhat.obj=cbhhat.obj))
}

summary_stat=function(Delta,time,X,beta,eta,xi,phi,psi,cbhhat.obj){
  ### input:  estimations at local site, make sure the dimension accords what is defined in the paper
  ### output: grad: gradient; 
  ###         Hessian: hessian
  Q = length(cbhhat.obj$weights)
  n = ncol(X)
  d_x = nrow(X)
  d_psi = length(eta)
  d_phi = nrow(phi)
  Z = as.numeric(t(xi)%*%phi+t(beta)%*%X)
  grad.lambda = matrix(0,d_psi,length(time))
  for (d in 1:d_psi){
    grad.lambda[d,] = (time/2)*(cbhhat.obj$weights%*%(cbhhat.obj$bh.eva*cbhhat.obj$psi.eva[,d,]))#colSums(cbhhat.obj$psi.eva[,d,]*(cbhhat.obj$bh.eva*cbhhat.obj$weights*time/2))
  }
  grad1 = -rowSums(sweep(psi,MARGIN=2,Delta,'*')-sweep(grad.lambda,exp(Z),MARGIN=2,'*'))
  grad2 = -rowSums(sweep(phi,MARGIN=2,Delta,'*')-sweep(phi,exp(Z)*cbhhat.obj$int,MARGIN=2,'*'))
  grad3 = -rowSums(sweep(X,MARGIN=2,Delta,'*')-sweep(X,exp(Z)*cbhhat.obj$int,MARGIN=2,'*'))
  grad = c(grad1,grad2,grad3)
  H1.1_ = array(0,dim=c(d_psi,d_psi,Q))
  H1.2_ = array(0,dim=c(d_psi,d_phi,Q))
  H1.3_ = array(0,dim=c(d_psi,d_x,Q))
  for (q in 1:Q){
    H1.1_[,,q] = cbhhat.obj$psi.eva[q,,]%*%
      diag(as.vector(time/2*exp(Z)*cbhhat.obj$weights[q]*exp(t(eta)%*%cbhhat.obj$psi.eva[q,,])),nrow=n,ncol=n)%*%
      t(cbhhat.obj$psi.eva[q,,])
    H1.2_[,,q]=cbhhat.obj$psi.eva[q,,]%*%
      diag(as.vector(time/2*exp(Z)*cbhhat.obj$weights[q]*exp(t(eta)%*%cbhhat.obj$psi.eva[q,,])),nrow=n,ncol=n)%*%
      t(phi)
    H1.3_[,,q]=cbhhat.obj$psi.eva[q,,]%*%
      diag(as.vector(time/2*exp(Z)*cbhhat.obj$weights[q]*exp(t(eta)%*%cbhhat.obj$psi.eva[q,,])),nrow=n,ncol=n)%*%
      t(X)
  }
  H1.1=rowSums(H1.1_,dims=2)
  H1.2=rowSums(H1.2_,dims=2)
  H1.3=rowSums(H1.3_,dims=2)
  H2.2=phi%*%
    diag(as.vector(exp(Z)*cbhhat.obj$int),nrow=n,ncol=n)%*%
    t(phi)
  H2.3=phi%*%
    diag(as.vector(exp(Z)*cbhhat.obj$int),nrow=n,ncol=n)%*%
    t(X)
  H3.3=X%*%
    diag(as.vector(exp(Z)*cbhhat.obj$int),nrow=n,ncol=n)%*%
    t(X)  
  H=rbind(cbind(H1.1,H1.2,H1.3),cbind(t(H1.2),H2.2,H2.3),cbind(t(H1.3),t(H2.3),H3.3))
  return(list(grad=1/n*grad,Hessian=1/n*H))
}
rmut = function(X,y){
  ### multiply each row of X with y, elementwise
  ### input X: a*b-dim matrix; y: b-dim array
  ### output a*b-dim matrix
  return(t(apply(X, 1, function(x){x*as.array(y)})))
}
cmut = function(X,y){
  ### multiply each column of X with y, elementwise
  ### input X: a*b-dim matrix; y: a-dim array
  ### output a*b-dim matrix
  return(apply(X, 2, function(x){x*y}))
}
get_true_bh = function(time,lambda,rho,dist_family){
  ### construct baseline hazard given parameters
  cbh=bh=length(time)
  if(dist_family=='weibull'){
    hazard = lambda*rho*time^(rho-1)
  }
  if(dist_family=='gompertz'){
    hazard = lambda*exp(rho*time)
  }
  return(data.frame(time=time,hazard=hazard))
}
com_psi0 = function(t,knot){
  if (length(t)==1){
    j = which(knot>=t)[1]-1
    output = rep(0,2*(length(knot)-1))
    output[2*(j)-1] = 1
    output[2*(j)] = t
  }else{
    output = matrix(0,length(t),2*(length(knot)-1))
    colnames(output)=paste0('psi',1:(2*(length(knot)-1)))
    for (j in 1:(length(knot)-1)){
      idx = t>=knot[j]&t<knot[j+1]
      output[idx,paste0('psi',((2*j)-1))] = 1
      output[,paste0('psi',(2*j))] = 0
      output[idx,paste0('psi',(2*j))] = t[idx]
    }
  }
  return(output)
}

com_psi = function(t,knot){
  output = cbind(1,t,matrix(0,length(t),length(knot)-2))
  colnames(output)=paste0('psi',1:ncol(output))
  for (j in 2:(length(knot)-1)){
    output[,paste0('psi',j+1)] = sapply(t, function(x){max(x-knot[j],0)})
  }
  return(output)
}

bhhat = function(t,eta,knot){
  psi = com_psi(t,knot)
  return(exp(psi%*%eta))
}
cbhhat = function(t,eta,knot,n.nodes=15){
  ### Approximate cumulative baseline hazard using gauss quadrature. 
  ### input: t: a vector of upper bounds of the integral of baseline hazard;
  ###        eta: coefficients corresponding to piece-wise exp approximation of baseline hazard
  ###        knot: knots used to for piece-wise exp approximation of baseline hazard
  ###        n.nodes: number of nodes used to approximate cumulative baseline hazard with gauss quadrature
  ### output: int: a vector of cumulative baseline hazard approximated by gauss quadrature with different upper bound specified by t
  ###         weights: weights used in gq (needed in computing grad and hessian)
  ###         bh.eva: the baseline hazard evaluated at nodes (needed in computing grad and hessian)
  ###         psi.eva: psi evaluted at knots
  myf = function(x){
    return(bhhat(x,eta=eta,knot=knot))
  }
  int = rep(0,length(t))
  weights = bh.eva = matrix(0,nrow=n.nodes,ncol=length(t))
  psi.eva = array(0,dim=c(n.nodes,length(eta),length(t)))
  for (i in 1:length(t)) {
    gq.res = gq.approx(myf,lower=0,upper=t[i],n.node=15)
    int[i] = gq.res$int
    bh.eva[,i]=gq.res$bh.eva
    psi.eva[,,i]=com_psi(gq.res$nodes,knot)
  }
  weights=gq.res$weights
  return(list(int=int,weights=weights,bh.eva=bh.eva,psi.eva=psi.eva))
}
gq.approx = function(myf,lower,upper,n.node=15){
  ### function performing gauss quadrature (gq) to approximate integral, 
  ### the idea of gq is change the integration into sum of weighted function value (evaluated at some nodes)
  ### see: https://coast.nd.edu/jjwteach/www/www/30125/pdfnotes/lecture16_19v17.pdf for example
  ### input: myf: the integral
  ###        lower: lower bound of integral
  ###        upper: upper bound of the integral
  ###        n.nodes: number of nodes taken in gq, only use the default value=15
  ### output: int: the value of the integral by approximation
  ###         weights: weights used in gq (needed in computing grad and hessian)
  ###         bh.eva: the baseline hazard evaluated at nodes (needed in computing grad and hessian)
  ###         nodes: the nodes used in gq (needed in computing grad and hessian)
  gaussquad = gauss.quad(n.node,kind='legendre')
  b = upper
  a = lower
  nodes = a+(b-a)/2*(gaussquad$nodes+1)
  eva= myf(nodes)
  integral = sum((b-a)/2*gaussquad$weights*eva)
  return(list(int=integral,weights=gaussquad$weight,bh.eva=eva,nodes=nodes))
}

surv_local = function(time,event,x,calendar_time, num_bh_knot=2, d_phi=3,
                      mean = NULL, sd = NULL){
  ### the main function for local parameter estimation and summary statistics derivation
  ### input: time: in numeric, survival time; 
  ###        x: in numeric, clinical feature (if value is 0 and included in missing_col, will construct missing indicator for the corresponding column) ; 
  ###        event: in binary, censorship (1 for not censored); 
  ###        calendar_time: in numeric, calendar_time, used to construct phi 
  ###        num_bh_knot: in int, number of nodes of piecewise-exponental function we use to approximate baseline hazard, num_bh_knot=dimension of psi, we can simply use the defalut value 2 nodes, corresponding to dim(psi)=4
  ###        d_phi: in tin, number of degrees of freedom of cubic spline for calendar time, we can use default value, 3
  ### output: a list containing estimations for all parameters and summary statistics, including gradient and Hessian
  # add cubic baseline for calendar time
  phi = data.frame(ns(calendar_time,df=d_phi))
  colnames(phi) = paste0('phi',1:ncol(phi))
  
  "
  if (is.null(mean)){
    for (j in 1:ncol(x)) {
      x[,j] <- (x[,j] - mean(x[,j])) 
    }
  }else{
    for (j in 1:ncol(x)) {
      if (!is.na(mean[j])){
        x[,j] <- x[,j] - mean[j]
      }
    }
  }
  if (is.null(sd)){
    for (j in 1:ncol(x)) {
      x[,j] <- x[,j] / sd(x[,j])
    }
  }else{
    for (j in 1:ncol(x)) {
      if (!is.na(sd[j])){
        x[,j] <- x[,j] / sd[j]
      }
    }
  }
  
  for (j in 1:d_phi) {
    phi[,j] <- (phi[,j] - mean(phi[,j])) / sd(phi[,j])
  }
  "
  
  # interaction with covariate
  X_ob0 = x
  X_ob = X_ob0
  for (j in 1:d_phi) {
    inter_basis = phi[,j]*X_ob0
    colnames(inter_basis) = paste0(colnames(phi)[j],".",colnames(X_ob0))
    X_ob = X_ob %>% mutate(inter_basis)
  }
  # interaction with missing indicator
  dat.obs = mutate(cbind(X_ob, phi),y=time,failed=event,calendar_time=calendar_time)
  local_est = get_local_est(dat.obs,num_bh_knot,d_phi)
  ### use estimations to compute summary statistics
  sum_stat = summary_stat(Delta=dat.obs$failed,
                          time=dat.obs$y,
                          X=local_est$X,
                          beta=local_est$beta,eta=local_est$eta,
                          xi=local_est$xi,phi=local_est$phi,
                          psi=local_est$psi,cbhhat.obj=local_est$cbhhat.obj)
  return(c(local_est,sum_stat))
}

