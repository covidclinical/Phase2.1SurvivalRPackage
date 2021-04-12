FUN.predict.cox <- function(newZ,beta,X,d,Z, t0, do_RD=T)
{
  # get rid of missing data
  cc <- complete.cases(cbind(X,d,Z));
  X <- X[cc];
  d <- d[cc];
  Z <- as.matrix(as.matrix(Z)[cc,]);
  
  # estimate cumulative baseline hazard
  dtimes <- sort(X[d==1]);
  # rows = death times, cols = subjs, 1 if at risk, 0 else
  atrisk <- outer(unique(dtimes),X,FUN="<=");
  
  scores <- as.matrix(Z)%*%as.matrix(beta);
  htemp <- 1/apply(atrisk,1,function(row)
  { sum(row*exp(scores)); });
  tied <- dtimes[duplicated(dtimes)];
  #for(tie in tied):  Coding error found with analysis BEST data
  for(tie in unique(tied))
  {
    # breslow method for dealing with ties
    # from p. 116 of kalbfleisch and prentice 2002
    htemp[unique(dtimes)==tie] <- htemp[unique(dtimes)==tie]*
      sum(dtimes==tie);
  }
  hazard <- cumsum(htemp);
  
  if (do_RD)
  {
    Lambda0 <- approx(unique(dtimes),hazard,t0)$y;
    # estimated predicted probabilities
    predicted <- exp(-exp(as.matrix(newZ)%*%as.matrix(beta))*
                       Lambda0);
  }else{
    # IRD
    time1 <- c(0,unique(dtimes))
    time2 <- c(unique(dtimes),1e+8)
    
    aa <- pmin(time2,tL)-pmin(time1,tL)
    bb <- pmin(time2,tU)-pmin(time1,tU)
    
    yy <- exp( -VTM(c(0,hazard),nrow(newZ))*as.vector(exp(as.matrix(newZ)%*%as.matrix(beta))) )
    predicted <- colSums(t(yy)*(bb-aa))/(tU-tL)
    # IRD
  }
  
  return(predicted);
}