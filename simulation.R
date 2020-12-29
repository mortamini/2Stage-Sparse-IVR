#
rm(list=ls())
#
setwd(choose.dir())
library(lars)
library(MASS)
library(grpreg)
source("EPIV2S.R")
source("LassoIV2S.R")
source("SCADIV2S.R")
source("eval-model.R")
source("gen-IV.R")
source("cv-model.R")
source("progress.R")
source("tcinv.R")
#
p=100
q=200
n=30
beta=c(rep(1,7),rep(0,85),rep(-0.5,8))
Gam=matrix(c(rep(0.01,500),rep(0,18500),rep(-0.005,1000)),q,p)
Z=matrix(rbinom(n*q,1,0.5),n,q)
#
pmat=combn(seq(0.1,0.9,0.2),2)
BB=100
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
	model1 <- EPIV2S(X,y,Z,trace=0,post=TRUE,
			intercept=c(FALSE,FALSE),init=list(model=model2))
	cat("evaluation ... \n")
	ev[[1]][bb,1:5] = eval.model(model1)
	cat("cross-validation ... \n")
	ev[[1]][bb,6] = cv.model("EPIV2S",trace=0,post=TRUE,
			intercept=c(FALSE,FALSE),init=list(model=model2))
 }, error=function(e){cat("Error in 2S EP:",
		conditionMessage(e), "\n")})

}
#-------------------------------------------------------
#
RE1=ev[[1]]/ev[[2]]
RE2=ev[[1]]/ev[[3]]
#
RE1[is.na(RE1) | is.nan(RE1)]=1
RE2[is.na(RE2) | is.nan(RE2)]=1
#
for(j in 1:6){
	RE1[!is.finite(RE1[,j])| RE1[,j]>quantile(RE1[,j],0.9)+1e-3 |RE1[,j]<quantile(RE1[,j],0.1)-1e-3,j]=NA
	RE2[!is.finite(RE2[,j]) | RE2[,j]>quantile(RE2[,j],0.9)+1e-3 |RE2[,j]<quantile(RE2[,j],0.1)-1e-3,j]=NA
}
#
#
par(mfrow=c(2,3))
boxplot(list("2S.EP/2S.LASSO"=RE1[,1],
"2S.EP/2S.SCAD"=RE2[,1]),
main=expression(paste("Relative FPR for estimator of ",beta)),
outline=FALSE)

boxplot(list("2S.EP/2S.LASSO"=RE1[,2],
"2S.EP/2S.SCAD"=RE2[,2]),
main=expression(paste("Relative FPR for estimator of ",Gamma)),
outline=FALSE)

boxplot(list("2S.EP/2S.LASSO"=RE1[,3],
"2S.EP/2S.SCAD"=RE2[,3]),
main=expression(paste("Relative FNR for estimator of ",beta)),
outline=FALSE)

boxplot(list("2S.EP/2S.LASSO"=RE1[,4],
"2S.EP/2S.SCAD"=RE2[,4]),
main=expression(paste("Relative FNR for estimator of ",Gamma)),
outline=FALSE)

boxplot(list("2S.EP/2S.LASSO"=RE1[,5],
"2S.EP/2S.SCAD"=RE2[,5]),
main="Relative RMSPE",outline=FALSE)

boxplot(list("2S.EP/2S.LASSO"=RE1[,6],
"2S.EP/2S.SCAD"=RE2[,6]),
main="Relative 3-fold CV",outline=FALSE)
#-------------------------------------------------------
