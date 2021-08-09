library(glmnet)
library(CVXR)

##########################################################################
###### The following function computes the lasso estimator
# Compute the Lasso estimator:
# - If lambda is given, use glmnet and standard Lasso
# - If lambda is set to the character string "CV", then glmnet with
#   lambda selected by cross-validation is used with 'lambda.1se'
# - If lambda is set to the character string "CV.min", then glmnet with
#   lambda selected by cross-validation is used with 'lambda.min'
##########################################################################
LASSO <- function(X, y, lambda=NULL, intercept=TRUE){
  p <- ncol(X)
  
  htheta <- if (lambda == 'CV') {
    outLas <- cv.glmnet(X, y, family='binomial', alpha=1, intercept=intercept)
    as.vector(coef(outLas, s=outLas$lambda.1se))
  } else if (lambda == 'CV.min'){
    outLas <- cv.glmnet(X, y, family='binomial', alpha=1, intercept=intercept)
    as.vector(coef(outLas, s=outLas$lambda.min))
  } else {
    outLas <- cv.glmnet(X, y, family='binomial', alpha=1, intercept=intercept)
    as.vector(coef(outLas, s=lambda))
  }
  
  if (intercept == TRUE) {
    return(htheta)
  } else {
    return(htheta[2:(p+1)])
  }
}

##########################################################################
###### Generate beta_star
# Arguments
# - B: coefficients of L groups, dim (p, L)
# - Sigma: covariates matrix, dim (p, p)
# - delta: default 0; the ridge penalty parameter; 
### if B^{T}\Sigma B is nearly singular, we may set delta to be non-zero, like 0.5.
##########################################################################
beta_star <- function(B, Sigma, delta=0){
  L <- dim(B)[2]
  Gamma <- t(B)%*%Sigma%*%B
  opt.weight <- rep(NA, L)
  # Problem Definition
  v <- Variable(L)
  Diag.matrix <- diag(eigen(Gamma)$values)
  for (ind in 1:L){
    Diag.matrix[ind, ind] <- max(Diag.matrix[ind, ind], 0.001)
  }
  Gamma.positive <- eigen(Gamma)$vectors %*% Diag.matrix %*% t(eigen(Gamma)$vectors)
  objective <- Minimize(quad_form(v, Gamma.positive + diag(delta, L)))
  constraints <- list(v>=0, sum(v)==1)
  prob.weight <- Problem(objective, constraints)
  if (is_dcp(prob.weight)){
    # Problem Solution
    result <- solve(prob.weight)
    opt.status <- result$status
    opt.sol <- result$getValue(v)
  }
  for (l in 1:L){
    opt.weight[l] <- opt.sol[l]*(abs(opt.sol[l])>10^{-8})
  }
  return(list(beta.est=B%*%opt.weight,weight=opt.weight))
}


