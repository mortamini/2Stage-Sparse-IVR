#
rm(list=ls())
#
setwd(choose.dir())
library(lars)
library(MASS)
library(hdm)
library(grpreg)
source("TSEPIV.R")
source("LassoIV2S.R")
source("SCADIV2S.R")
source("rlassoIV2.R")
#source("LassoIVSW.R")
source("eval-model.R")
source("gen-IV.R")
source("cv-model.R")
#
p=100
q=200
n=30
beta=c(rep(1,7),rep(0,85),rep(-0.5,8))
Gam=matrix(c(rep(0.01,500),rep(0,18500),rep(-0.005,1000)),q,p)
Z=matrix(rbinom(n*q,1,0.5),n,q)
#
BB=1000
#
ev=list()
ev[[1]]=ev[[2]]=ev[[3]]=ev[[4]]=matrix(nrow=BB,ncol=6)
ev[[5]]=ev[[6]]=ev[[7]]=ev[[8]]=matrix(nrow=BB,ncol=6)
for(bb in 1:BB){
	gen<-gen.IV(beta,Gam,Z,s02=0.5,tau02=0.1,
		int1=0.1,int2=1)
	y=gen$y
	X=gen$X
	#
	#---------------------------------------------------
tryCatch({
	cat(paste("iteration",bb,"2S-LASSO starts \n"))
	model2=LassoIV2S(X,y,Z)
	cat("evaluation ... \n")
	ev[[2]][bb,1:5] = eval.model(model2)
	cat("cross-validation ... \n")
	ev[[2]][bb,6] = cv.model("LassoIV2S")
  }, error=function(e){cat("Error in 2S Lasso:",
		conditionMessage(e), "\n")})
#	#---------------------------------------------------
tryCatch({
	cat(paste("iteration",bb,"2S-SCAD starts \n"))
	model3=SCADIV2S(X,y,Z)
	cat("evaluation ... \n")
	ev[[3]][bb,1:5] = eval.model(model3)
	cat("cross-validation ... \n")
	ev[[3]][bb,6] = cv.model("SCADIV2S")
  }, error=function(e){cat("Error in 2S SCAD:",
		conditionMessage(e), "\n")})
	#---------------------------------------------------
	#---------------------------------------------------
	cat(paste("iteration",bb,"2S-EP starts \n"))
tryCatch({
	model1 <- EPIV2S(X,y,Z,trace=1,post=FALSE,
			intercept=TRUE,init=model2)
	cat("evaluation ... \n")
	ev[[1]][bb,1:5] = eval.model(model1)
	cat("cross-validation ... \n")
	ev[[1]][bb,6] = cv.model("TSEPIV",post=FALSE,
			intercept=TRUE,init=LassoIV2S,
			init.criteria="AICc")
  }, error=function(e){cat("Error in 2S EP:",
		conditionMessage(e), "\n")})

}
#-------------------------------------------------------
#
RE1=ev[[2]]/ev[[1]]
RE2=ev[[3]]/ev[[1]]
#
RE1[is.na(RE1) | is.nan(RE1)]=1
RE1[!is.finite(RE1)| RE1>4 |RE1<0.25]=NA
RE2[is.na(RE2) | is.nan(RE2)]=1
RE2[!is.finite(RE2) | RE2>4 |RE2<0.25]=NA
#
#
par(mfrow=c(2,3))
boxplot(list("2S.LASSO/2S.EP"=RE1[,1],
"2S.SCAD/2S.EP"=RE2[,1]),
main=expression(paste("Relative FPR for estimator of ",beta)),
outline=FALSE)

boxplot(list("2S.LASSO/2S.EP"=RE1[,2],
"2S.SCAD/2S.EP"=RE2[,2]),
main=expression(paste("Relative FPR for estimator of ",Gamma)),
outline=FALSE)

boxplot(list("2S.LASSO/2S.EP"=RE1[,3],
"2S.SCAD/2S.EP"=RE2[,3]),
main=expression(paste("Relative FNR for estimator of ",beta)),
outline=FALSE)

boxplot(list("2S.LASSO/2S.EP"=RE1[,4],
"2S.SCAD/2S.EP"=RE2[,4]),
main=expression(paste("Relative FNR for estimator of ",Gamma)),
outline=FALSE)

boxplot(list("2S.LASSO/2S.EP"=RE1[,5],
"2S.SCAD/2S.EP"=RE2[,5]),
main="Relative RMSPE",outline=FALSE)

boxplot(list("2S.LASSO/2S.EP"=RE1[,6],
"2S.SCAD/2S.EP"=RE2[,6]),
main="Relative 3-fold CV",outline=FALSE)
#-------------------------------------------------------
