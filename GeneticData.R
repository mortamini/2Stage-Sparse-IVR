#
rm(list=ls())
#
setwd(choose.dir())
load(file="geneticdata.Rdata")
#
#---------------------------------------------------
p=dim(X)[2]
q=dim(Z)[2]
n=length(y)
#
library(lars)
library(MASS)
library(hdm)
library(grpreg)
source("EPIV2S.R")
source("LassoIV2S.R")
source("SCADIV2S.R")
source("eval-model.R")
source("gen-IV.R")
source("cv-model.R")
source("progress.R")
#---------------------------------------------------
	model2=LassoIV2S(X,y,Z,method="grpreg",trace=2)
	cv2 = cv.model("LassoIV2S",method="grpreg",trace=2)
	model3=SCADIV2S(X,y,Z)
	cv3 = cv.model("SCADIV2S")
	model1 <- EPIV2S(X,y,Z,trace=2,post=FALSE,
			intercept=TRUE,init=list(model2,p0=0.5,pi0=0.5),
			eps1=1e-4,eps2=1e-4)
	cv1 = cv.model("TSEPIV",trace=2,post=FALSE,
			intercept=TRUE,init=list(model2,p0=0.5,pi0=0.5),
			eps1=1e-4,eps2=1e-4)
#---------------------------------------------------
#---------------------------------------------------
#--------------------LassoIV2S------------------------
#---------------------------------------------------
betahat = model2$coefficients
Gammahat = model2$gamma
Xhat=cbind(1,Z)%*%Gammahat
yhat=cbind(1,Xhat)%*%betahat
R2=1-sum((yhat-y)^2)/sum((y-mean(y))^2)
R2
BIC=n*log(sum((yhat-y)^2)/n) + p*log(n)
BIC
#---------------------------------------------------
#---------------------------------------------------
#--------------------SCADIV2S------------------------
#---------------------------------------------------
betahat = model3$coefficients
Gammahat = model3$gamma
Xhat=cbind(1,Z)%*%Gammahat
yhat=cbind(1,Xhat)%*%betahat
R2=1-sum((yhat-y)^2)/sum((y-mean(y))^2)
R2
BIC=n*log(sum((yhat-y)^2)/n) + p*log(n)
BIC
#---------------------------------------------------
#---------------------------------------------------
#--------------------TSEPIV------------------------
#---------------------------------------------------
betahat = model1$coefficients
Gammahat = model1$gamma
Xhat=cbind(1,Z)%*%Gammahat
Xhat=Z%*%Gammahat
yhat=cbind(1,Xhat)%*%betahat
yhat=Xhat%*%betahat
R2=1-sum((yhat-y)^2)/sum((y-mean(y))^2)
R2
BIC=n*log(sum((yhat-y)^2)/n) + p*log(n)
BIC
#---------------------------------------------------
sparseplot<-function(y,x,...){
plot(y~x,type="n",...)
n=length(x)
for(i in 1:n){
lines(c(x[i],x[i]),c(y[i],0))
}
}
resize.win <- function(Width=6, Height=6)
{
        # works for windows
    dev.off(); # dev.new(width=6, height=6)
    windows(record=TRUE, width=Width, height=Height)
}
#---------------------------------------------------
tt=1:(p+1)
ss=1:(q+1)
resize.win(40,20)
par(mfrow=c(1,2))
qgamma=-1*(Gammahat==0)+1*(Gammahat!=0)
sparseplot(betahat,tt,xlab="Gene number",ylab="Linear effect")
image(ss,tt,qgamma,col=c("white","gray","black"),
xlab="SNP number",ylab="Gene number")
#----------------