Est.ALASSO.GLMNET=function (data, BIC.factor = 0.1, fam0 = "binomial", w.b = NULL, Wi, lambda.grid = 10^seq(-4,0,0.01)) 
{  if(is.null(Wi)==1){Wi=rep(1,dim(data)[1])}
  if (fam0 != "Cox") {
    data = as.matrix(data)
    y = data[, 1]
    x = data[, -1, drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
    if (is.null(w.b)) {
      bini = as.vector(coef(glmnet(x,y, family=fam0, weights=Wi,alpha=0,standardize=F,lambda=0.1)))
      w.b = 1/abs(bini[-1])
      gc()
    }
    tmpfit = glmnet(x = x, y = y, family = fam0, penalty.factor = w.b, weights=Wi,
                    alpha = 1, lambda = lambda.grid, intercept = T)
    gc()
    N.adj = dim(x)[1]
  }
  else {
    data = as.matrix(data)
    y = cbind(time = data[, 1], status = data[, 2])
    x = data[, -c(1, 2), drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
    if (is.null(w.b)) {
      bini = as.vector(coef(glmnet(x, y, family = "cox",  weights=Wi,
                      alpha = 0, lambda = 0.1)))
      w.b = 1/abs(bini)
      gc()
    }
    tmpfit = glmnet(x, y, family = "cox", penalty.factor = w.b, weights=Wi,
                    alpha = 1, lambda = lambda.grid)
    gc()
    N.adj = sum(y[, 2])
  }
  dev = deviance(tmpfit)
  BIC.lam = dev + min(N.adj^0.1, log(N.adj)) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.modBIC = c(tmpfit$a0[m.opt], tmpfit$beta[, m.opt])
  lamhat.modBIC = tmpfit$lambda[m.opt]
  BIC.lam = dev + log(N.adj) * tmpfit$df
  m.opt = which.min(BIC.lam)
  bhat.BIC = c(tmpfit$a0[m.opt], tmpfit$beta[, m.opt])
  lamhat.BIC = tmpfit$lambda[m.opt]
  
  
  tLL <- tmpfit$nulldev - deviance(tmpfit)
  k <- tmpfit$df
  n <- tmpfit$nobs
  AIC.lam <- -tLL+2*k+2*k*(k+1)/(n-k-1)
  
  m.opt=which.min(AIC.lam)
  bhat.AIC=c(tmpfit$a0[m.opt], tmpfit$beta[,m.opt])
  lamhat.AIC = tmpfit$lambda[m.opt]
  
  return(list(bhat.BIC = bhat.BIC, bhat.modBIC = bhat.modBIC, bhat.AIC=bhat.AIC,
              lambda.BIC = lamhat.BIC, lambda.modBIC = lamhat.modBIC, lambda.AIC=lamhat.AIC))
}