MySum <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
{
  if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
  if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
  pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
  
  if(is.null(Vi)){return(pos)}else{
    Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
    out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
    out[pos!=0,] <- Vi[pos,]
    if(is.null(dim(Vi))) out <- c(out)
    return(out) ## n.y x p
  }
}

PI.k.FUN <- function(tt,ebzi,xi,zi,k0=0,vi=NULL)
{
  out = ebzi; pz = ncol(zi); nn=length(out)
  if(k0==1){out=out*zi}
  if(k0==2){out = out*zi[,rep(1:pz,pz)]*zi[,rep(1:pz,rep(pz,pz))]}
  as.matrix(MySum(tt,"<=",xi,Vi=out)/nn)
}

Score.A.FUN <- function(data,betahat,rtn="score")
{
  xi = data[,1]; di = data[,2]; zi = data.matrix(data[,-(1:2),drop=F]); ebzi = c(exp(zi%*%betahat)); pz = ncol(zi); nn = length(xi);
  tmpind = di==1; tj = xi[tmpind]; zj = zi[tmpind,,drop=F]
  pi0.tj   = c(PI.k.FUN(tj,ebzi,xi,zi,k0=0)); pi1.tj   = PI.k.FUN(tj,ebzi,xi,zi,k0=1)
  smat = zj - pi1.tj/pi0.tj
  score = apply(smat,2,sum)/dim(data)[1]
  if(rtn=="score"){
    return(score)
  }else if (rtn=='Score+A'){
    Ahat = PI.k.FUN(tj,ebzi,xi,zi,k0=2)*pi0.tj - pi1.tj[,rep(1:pz,pz)]*pi1.tj[,rep(1:pz,rep(pz,pz))]
    Ahat = -matrix(apply(Ahat/pi0.tj^2,2,sum),ncol=pz)/dim(data)[1]
    return(list("score"=score,"neg_info_mat"=Ahat))
  }else if (rtn=='Score+A+Approx'){
    Ahat = -t(smat) %*% smat/dim(data)[1]
    return(list("score"=score,"neg_info_mat"=Ahat))
  }
}

iteration.fun.cox = function(dat.list,bini,kk.list,rtn='Score+A+Approx'){
  K = length(dat.list)
  ScoreA.list = lapply(1:K,function(kk){Score.A.FUN(dat.list[[kk]],bini,rtn)})
  Uini.list = sapply(kk.list,function(kk){ScoreA.list[[kk]]$score})
  Aini.list = lapply(1:K,function(kk){ScoreA.list[[kk]]$neg_info_mat});Ahat.ini = Reduce("+",Aini.list)/K
  bhat.list = -solve(Ahat.ini)%*%Uini.list + bini;
  list("b.k"=bhat.list,"Ahat"=Ahat.ini)
}
