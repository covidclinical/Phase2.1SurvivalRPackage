Est.ALASSO.GLMRIDGE=function (data, BIC.factor = 0.1, fam0 = "binomial",  Wi, lambda.grid = 10^seq(-4,0,0.01)) 
{  if(is.null(Wi)==1){Wi=rep(1,dim(data)[1])}
  if (fam0 != "Cox") {
    data = as.matrix(data)
    y = data[, 1]
    x = data[, -1, drop = F]
    nn = length(y)
    pp = ncol(x)

      bini = as.vector(coef(glmnet(x,y, family=fam0, weights=Wi,alpha=0,standardize=F,lambda=1/nn)))

  }
  else {
    data = as.matrix(data)
    y = cbind(time = data[, 1], status = data[, 2])
    x = data[, -c(1, 2), drop = F]
    nn = length(y)
    pp = ncol(x)
    gc()
      
      bini = as.vector(coef(glmnet(x, y, family = "cox",  weights=Wi,
                      alpha = 0, lambda = 1/nn)))

}
  
  return(list(bhat=bini))
}