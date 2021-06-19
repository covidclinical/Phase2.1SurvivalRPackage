
X = iris[,c(1:4)]
n=dim(X)[1]
Y=data.matrix(cbind(1,X))%*%c(1,0.5,0.2,0.1,0.1)+rnorm(n, 0,1)
junk1=lm(Y~data.matrix(X))
coef(junk1)[-1]

mu = colMeans(X)
Xpca = prcomp(X, center=F, scale=F)
Xpca$rotation
PC0=data.matrix(X)%*%Xpca$rotation

Xcov <- cov(X)
Xeigen <- eigen(Xcov)
pca.loadings <- Xeigen$vectors
PC=as.matrix(X)%*%pca.loadings
cov(PC)

junk2=lm(Y~Xpca$x)
coef(junk2)[-1]%*%t(Xpca$rotation)


junk3=lm(Y~PC)
coef(junk3)[-1]%*%t(pca.loadings)

