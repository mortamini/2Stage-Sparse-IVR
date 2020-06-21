LassoIV2S<-function(X,y,Z,trace=0,intercept=TRUE,criteria="AICc",method="lars"){
	#---------------------------------------------------
	if(nrow(X)!=length(y) | nrow(X)!=nrow(Z)) stop("number of samples 
			of X, y and Z are not the same!")
	#---------------------------------------------------	
	n=nrow(X)
	p=ncol(X)
	q=ncol(Z)
	Z1=Z
	if(method=="grpreg") intercept=TRUE
	#---------------------------------------------------
	if(trace>=1) cat("Stage I ...\n")
	Gh=matrix(nrow=q+intercept,ncol=p)
	for(j in 1:p){
		if(trace==2) cat("Covariate model",j,"...\n")
		if(method=="lars"){
		tryCatch({
		if(q>500){
			model2=lars(Z,X[,j],type = "lasso",intercept =intercept,use.Gram=FALSE)
		}else{
			model2=lars(Z,X[,j],type = "lasso",intercept =intercept)			
		}
		}, error=function(e){
			stop("An error occured in LARS ! Use method=grpreg instead. \n")
		})
		}else if(method =="grpreg"){
			model2=grpreg(Z, X[,j], nlambda=min(q,100), alpha=1, intercept=FALSE, 
					group=1:q, penalty="grLasso", family="gaussian",
					group.multiplier=rep(1,q))
		}else{ stop("method must be \"lars\" or \"grpreg\"")}
		gamhat=coef(model2)
		if(method=="lars"){
			g0=predict(model2,data.frame(t(rep(0,q))))$fit
			gamhat=cbind(g0,gamhat)
		}
		if(method =="grpreg") gamhat=t(gamhat)
		gamhat=gamhat[-1,]
		if(intercept) Z1=cbind(1,Z)
		Xhats=Z1%*%t(gamhat)
		resids=X[,j]-Xhats
		RSSs=apply(resids,2,function(x){sum(x^2)})
		dfs=apply(gamhat[,-1],1,function(x){sum(x != 0)})
		if(criteria=="AICc")
			CRTs=2*dfs+n*log(RSSs/n)+(2*dfs^2+2*dfs)/(n-dfs-1)
		if(criteria=="AIC")
			CRTs=2*dfs+n*log(RSSs/n)
		if(criteria=="BIC")
			CRTs=log(n)*dfs+n*log(RSSs/n)
		if(criteria=="BICc")
			CRTs=log(q)*dfs+n*log(RSSs/n)
		path = which.min(CRTs)
		Gh[,j]=gamhat[path,]
	}
	if(trace>=1) cat("Stage II ...\n")
	Xhat=Z1%*%Gh
	if(method=="lars"){
	tryCatch({
	if(p>500){
		model2=lars(Xhat,y,type = "lasso",intercept =intercept,use.Gram=FALSE)
	}else{
		model2=lars(Xhat,y,type = "lasso",intercept =intercept)
	}
	}, error=function(e){
			stop("An error occured in LARS ! Use method=grpreg instead. \n")
	})
	}else if(method=="grpreg"){
		model2=grpreg(Xhat,y, nlambda=min(p,100), alpha=1, intercept=FALSE, 
				group=1:p, penalty="grLasso", family="gaussian",
				group.multiplier=rep(1,p))
	}
	bethat=coef(model2)
	if(method=="lars"){
		b0=predict(model2,data.frame(t(rep(0,p))))$fit
		bethat=c(b0,bethat)
	}
	if(method=="grpreg") bethat=t(bethat)
	bethat=bethat[-1,]
	if(intercept) Xhat=cbind(1,Xhat)
	yhats=Xhat%*%t(bethat)
	resids=y-yhats
	RSSs=apply(resids,2,function(x){sum(x^2)})
	dfs=apply(bethat[,-1],1,function(x){sum(x != 0)})
		if(criteria=="AICc")
			CRTs=2*dfs+n*log(RSSs/n)+(2*dfs^2+2*dfs)/(n-dfs-1)
		if(criteria=="AIC")
			CRTs=2*dfs+n*log(RSSs/n)
		if(criteria=="BIC")
			CRTs=log(n)*dfs+n*log(RSSs/n)
		if(criteria=="BICc")
			CRTs=log(p)*dfs+n*log(RSSs/n)
	path = which.min(CRTs)
	bh=bethat[path,]
	sel = Gh!=0
	out=list(coefficients=bh,selected=sel,gamma=Gh)
	return(out)
}