mysurvfit.multi.fun=function(betahat, data, newdata, nm.x){
  tab <- data.frame(table(data[data$status == 1, "time"])) 
  y <- as.numeric(levels(tab[, 1]))[tab[, 1]] 
  d <- tab[, 2]       
  h0 <- rep(NA, length(y))
  for(l in 1:length(y))
  {
    h0[l] <- d[l] / sum(exp(data.matrix(data[data$time >= y[l],nm.x]) %*% betahat))
  }
  surv=exp(matrix(-c(exp(data.matrix(newdata[,nm.x]) %*% betahat)))%*%cumsum(h0))
  colnames(surv)=paste0("t", 1:dim(surv)[2])
  surv
}