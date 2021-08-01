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
source("MH.R")
#
p=300
q=400
n=50
beta=c(rep(1,7),rep(0,285),rep(-0.5,8))
Gam=matrix(c(rep(0.01,300),rep(0,118800),rep(-0.005,900)),q,p)
Z=matrix(rbinom(n*q,1,rbeta(n*q,3,7)),n,q)
#
pmat=combn(seq(0.1,0.9,0.2),2)
BB=1000
#
ev=list()
ev[[1]]=ev[[2]]=ev[[3]]=ev[[4]]=matrix(nrow=BB,ncol=7)
#ev[[5]]=ev[[6]]=ev[[7]]=ev[[8]]=matrix(nrow=BB,ncol=6)
for(bb in 1:BB){
	gen<-gen.IV(beta,Gam,Z,s02=0.5,tau02=0.1,
		int1=0.1,int2=1)
	y=gen$y
	X=gen$X
	#
	#---------------------------------------------------
tryCatch({
	start_time <- Sys.time()
	cat(paste("iteration",bb,"2S-LASSO starts \n"))
	model2=LassoIV2S(X,y,Z,method ="grpreg")
	cat("evaluation ... \n")
	ev[[2]][bb,1:5] = eval.model(model2)
	cat("cross-validation ... \n")
	ev[[2]][bb,6] = cv.model("LassoIV2S",method ="grpreg")
	end_time <- Sys.time()
	ev[[2]][bb,7] = end_time - start_time
 }, error=function(e){cat("Error in 2S Lasso:",
		conditionMessage(e), "\n")})
#	#---------------------------------------------------
tryCatch({
	start_time <- Sys.time()
	cat(paste("iteration",bb,"2S-SCAD starts \n"))
	model3=SCADIV2S(X,y,Z)
	cat("evaluation ... \n")
	ev[[3]][bb,1:5] = eval.model(model3)
	cat("cross-validation ... \n")
	ev[[3]][bb,6] = cv.model("SCADIV2S")
	end_time <- Sys.time()
	ev[[3]][bb,7] = end_time - start_time
  }, error=function(e){cat("Error in 2S SCAD:",
		conditionMessage(e), "\n")})
	#---------------------------------------------------
	#---------------------------------------------------
	cat(paste("iteration",bb,"2S-EP starts \n"))
tryCatch({
	start_time <- Sys.time()
	model1 <- EPIV2S(X,y,Z,trace=0,post=TRUE,
			intercept=c(FALSE,FALSE),init=list(model=model2))
	cat("evaluation ... \n")
	ev[[1]][bb,1:5] = eval.model(model1)
	cat("cross-validation ... \n")
	ev[[1]][bb,6] = cv.model("EPIV2S",trace=0,post=TRUE,
			intercept=c(FALSE,FALSE),init=list(model=model2))
	end_time <- Sys.time()
	ev[[1]][bb,7] = end_time - start_time
 }, error=function(e){cat("Error in 2S EP:",
		conditionMessage(e), "\n")})
}
	#---------------------------------------------------
	#---------------------------------------------------
#	cat(paste("iteration",bb,"MH starts \n"))
#tryCatch({
#	start_time <- Sys.time()
#	model4 <- AMH(X,y,Z,init=list(model=model2))
#	cat("evaluation ... \n")
#	ev[[4]][bb,1:5] = eval.model(model4)
#	cat("cross-validation ... \n")
#	ev[[4]][bb,6] = cv.model("AMH",init=list(model=model2))
#	end_time <- Sys.time()
#	ev[[4]][bb,7] = end_time - start_time
# }, error=function(e){cat("Error in MH:",
#		conditionMessage(e), "\n")})
#}
#-------------------------------------------------------
#
#RE1=ev[[1]]/ev[[2]]
#RE2=ev[[1]]/ev[[3]]
#
#RE1[is.na(RE1) | is.nan(RE1)]=1
#RE2[is.na(RE2) | is.nan(RE2)]=1
#
#for(j in 1:6){
#	RE1[!is.finite(RE1[,j])| RE1[,j]>quantile(RE1[,j],0.9)+1e-3 |RE1[,j]<quantile(RE1[,j],0.1)-1e-3,j]=NA
#	RE2[!is.finite(RE2[,j]) | RE2[,j]>quantile(RE2[,j],0.9)+1e-3 |RE2[,j]<quantile(RE2[,j],0.1)-1e-3,j]=NA
#}
#
#
par(mfrow=c(2,3))
boxplot(list("2S.EP"=ev[[1]][,1],
"2S.SCAD"=ev[[3]][,1],"2S.LASSO"=ev[[2]][,1]),
#,"M-H"=ev[[4]][,1]),
main=expression(paste("FPR for estimator of ",beta)),
outline=FALSE)

boxplot(list("2S.EP"=ev[[1]][,2],
"2S.SCAD"=ev[[3]][,2],"2S.LASSO"=ev[[2]][,2]),
#,"M-H"=ev[[4]][,2]),
main=expression(paste("FPR for estimator of ",Gamma)),
outline=FALSE)

boxplot(list("2S.EP"=ev[[1]][,3],
"2S.SCAD"=ev[[3]][,3],"2S.LASSO"=ev[[2]][,3]),
#,"M-H"=ev[[4]][,3]),
main=expression(paste("FNR for estimator of ",beta)),
outline=FALSE)

boxplot(list("2S.EP"=ev[[1]][,4],
"2S.SCAD"=ev[[3]][,4],"2S.LASSO"=ev[[2]][,4]),
#,"M-H"=ev[[4]][,4]),
main=expression(paste("FNR for estimator of ",Gamma)),
outline=FALSE)

#boxplot(list("2S.EP"=ev[[1]][,5],
#"2S.SCAD"=ev[[3]][,5],"2S.LASSO"=ev[[2]][,5]),
#,"M-H"=ev[[4]][,5]),
#main="RMSPE",outline=FALSE)

boxplot(list("2S.EP"=ev[[1]][,6],
"2S.SCAD"=ev[[3]][,6],"2S.LASSO"=ev[[2]][,6]),
#,"M-H"=ev[[4]][,6]),
main="3-fold CV",outline=FALSE)

boxplot(list("2S.EP"=ev[[1]][,7],
"2S.SCAD"=ev[[3]][,7],"2S.LASSO"=ev[[2]][,7]),
#,"M-H"=ev[[4]][,7]),
main="Computation time",outline=FALSE)

#-------------------------------------------------------
