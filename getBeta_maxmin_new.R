setwd("/Users/chuanhong/Documents/GitHub/Phase2.1SurvivalRPackage/FourCePhase2.1Survival/")

source("R/logistic_maximin.R")

# library(readr)
# dir="/Users/xuanwang/Desktop/output_FederatedLearning_06192021/"
# beta.maxmin=read_csv(paste0(dir,"beta.maxmin.table.toXin.csv"))
load("data/dat_beta_xTransX_sd.Rdata")
resnew=res
#load("~/Desktop/maximin/fromChuan/beta.maxmin.rda")
#load("~/Desktop/maximin/fromChuan/xTransX.Rdata")
#load("~/Desktop/maximin/fromChuan/betahat.port.deceased.rda")

# all the beta
covariates=names(resnew$VA3$beta)
beta.all=NULL
for (i in 1:18){
  betatemp=(resnew[[i]]$beta)[match(covariates, names(resnew[[i]]$beta) )]
  beta.all=cbind(beta.all, betatemp)
}
colnames(beta.all)=names(resnew)
rownames(beta.all)=covariates

# ind.norace=c(1,3,4,5,10)
# for (j in 1:length(ind.norace)){
#   beta.all[c('raceBlack','raceAsian','raceHispanic.and.Other'),ind.norace[j]]=NA
#   res[[ind.norace[j]]]$xx[,c('raceBlack','raceAsian','raceHispanic.and.Other')]=NA
#   res[[ind.norace[j]]]$xx[c('raceBlack','raceAsian','raceHispanic.and.Other'),]=NA
# }

index=which(apply(is.na(beta.all),2,sum)>0)
beta.complete=beta.all[,-index]
name.site=colnames(beta.complete)
name.var=row.names(beta.complete)

# all matrix, all sd
for (i in 1:18){
  ind=match(covariates,colnames(resnew[[i]]$xx))
  resnew[[i]]$xx=resnew[[i]]$xx[ind,ind]
  ind=match(covariates,names(resnew[[i]]$sd))
  resnew[[i]]$sd=resnew[[i]]$sd[ind]
}
res.complete=resnew

res.complete[[18]]=NULL
res.complete[[14]]=NULL
res.complete[[10]]=NULL
res.complete[[5]]=NULL
res.complete[[4]]=NULL
res.complete[[3]]=NULL
res.complete[[1]]=NULL

#### 
beta.est=matrix(0,nrow(beta.complete),ncol(beta.complete))
for (i in 1:ncol(beta.complete)){
# i=1
B = beta.complete[,-i]
meanx=sqrt(diag(res.complete[[i]]$xx / res.complete[[i]]$N.train)-(res.complete[[i]]$sd)^2)
Sigma = res.complete[[i]]$xx / res.complete[[i]]$N.train-matrix(meanx,ncol=1)%*%matrix(meanx,nrow=1)
# print(round(eigen(Sigma)$values,2))

# output=beta_star(B, Sigma, delta=0)
# reward0=t(output$beta.est)%*%Sigma%*%output$beta.est
# reward0=min( 2*t(B)%*%Sigma%*%output$beta.est-rep(t(output$beta.est)%*%Sigma%*%output$beta.est,ncol(B)) )
# output1=beta_star(B, Sigma, delta=0.5)
# reward=min( 2*t(B)%*%Sigma%*%output1$beta.est-rep(t(output1$beta.est)%*%Sigma%*%output1$beta.est,ncol(B)) )
# re1=((reward0-reward)/reward0)
# output2=beta_star(B, Sigma, delta=2)
# reward=min( 2*t(B)%*%Sigma%*%output2$beta.est-rep(t(output2$beta.est)%*%Sigma%*%output2$beta.est,ncol(B)) )
# re2=((reward0-reward)/reward0)
# if (re2<=0.05){output=output2}else if(re1<=0.05){output=output1}

output=beta_star(B, Sigma, delta=2)
beta.est[,i]=output$beta.est
}
rownames(beta.est)=rownames(beta.complete)
colnames(beta.est)=colnames(beta.complete)

#### 
library(MASS)
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

# race
# 1,3,4,5,10
beta.hat=matrix(0, nrow(beta.all),5)
for (l in 1:5){
k=index[[l]]
betatemp=beta.all[,k]
ind.na=which(is.na(betatemp))

beta.impute=matrix(0,nrow(beta.complete)-length(ind.na),ncol(beta.complete))
for (i in 1:ncol(beta.complete)){
# i=1
temp=ginv(res.complete[[i]]$xx[-ind.na,-ind.na]) %*% ( (res.complete[[i]]$xx[,ind.na])[-ind.na,] )
beta.impute[,i]=beta.complete[-ind.na,i]+apply(temp*VTM(beta.complete[ind.na,i],nrow(temp)),1,sum)
}
rownames(beta.impute)=name.var[-ind.na]
colnames(beta.impute)=name.site

B = beta.impute
meanx=sqrt(diag(resnew[[k]]$xx[-ind.na,-ind.na] / resnew[[k]]$N.train)-(resnew[[k]]$sd[-ind.na])^2)
Sigma = resnew[[k]]$xx[-ind.na,-ind.na] / resnew[[k]]$N.train-matrix(meanx,ncol=1)%*%matrix(meanx,nrow=1)
# print(round(eigen(Sigma)$values,2))
output=beta_star(B, Sigma, delta=2)
beta.hat[,l]=(output$beta.est)[match(name.var,rownames(output$beta.est))]
}
rownames(beta.hat)=name.var
colnames(beta.hat)=names(resnew)[c(1,3,4,5,10)]

# 14, 18
k=14
betatemp=beta.all[,k]
ind.na=which(is.na(betatemp))
beta.impute=matrix(0,nrow(beta.complete)-length(ind.na),ncol(beta.complete))
for (i in 1:ncol(beta.complete)){
  # i=1
  temp=ginv(res.complete[[i]]$xx[-ind.na,-ind.na]) %*% ((res.complete[[i]]$xx[,ind.na])[-ind.na])
  beta.impute[,i]=beta.complete[-ind.na,i]+temp*beta.complete[ind.na,i]
}
rownames(beta.impute)=name.var[-ind.na]
colnames(beta.impute)=name.site

B = beta.impute
meanx=sqrt(diag(resnew[[k]]$xx[-ind.na,-ind.na] / resnew[[k]]$N.train)-(resnew[[k]]$sd[-ind.na])^2)
Sigma = resnew[[k]]$xx[-ind.na,-ind.na] / resnew[[k]]$N.train-matrix(meanx,ncol=1)%*%matrix(meanx,nrow=1)
# print(round(eigen(Sigma)$values,2))
output=beta_star(B, Sigma, delta=2)
beta.hat=cbind(beta.hat,(output$beta.est)[match(name.var,rownames(output$beta.est))])
colnames(beta.hat)=names(res)[index[1:6]]

k=18
betatemp=beta.all[,k]
ind.na=which(is.na(betatemp))
beta.impute=matrix(0,nrow(beta.complete)-length(ind.na),ncol(beta.complete))
for (i in 1:ncol(beta.complete)){
  # i=1
  temp=ginv(res.complete[[i]]$xx[-ind.na,-ind.na]) %*% ((res.complete[[i]]$xx[,ind.na])[-ind.na])
  beta.impute[,i]=beta.complete[-ind.na,i]+temp*beta.complete[ind.na,i]
}
rownames(beta.impute)=name.var[-ind.na]
colnames(beta.impute)=name.site

B = beta.impute
meanx=sqrt(diag(resnew[[k]]$xx[-ind.na,-ind.na] / resnew[[k]]$N.train)-(resnew[[k]]$sd[-ind.na])^2)
Sigma = resnew[[k]]$xx[-ind.na,-ind.na] / resnew[[k]]$N.train-matrix(meanx,ncol=1)%*%matrix(meanx,nrow=1)
# print(round(eigen(Sigma)$values,2))
output=beta_star(B, Sigma, delta=2)
beta.hat=cbind(beta.hat,(output$beta.est)[match(name.var,rownames(output$beta.est))])
colnames(beta.hat)=names(resnew)[index]


temp=data.frame(cbind(beta.est,beta.hat))
beta.maximin=temp[,match(names(resnew),colnames(temp))]
meanbetaEU=apply(beta.maximin[,c(1,3,4,5,10)],1,mean,na.rm=T)
meanbetaUS=apply(beta.maximin[,-c(1,3,4,5,10)],1,mean,na.rm=T)
meanbetaall=apply(beta.maximin,1,mean,na.rm=T)
out=list('beta.site'=beta.all,'beta.maximin'=beta.maximin,'meanbetaall'=meanbetaall,'meanbetaEU'=meanbetaEU,'meanbetaUS'=meanbetaUS)
save(out,file='data/betamaximin.rda')

