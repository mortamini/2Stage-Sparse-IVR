SCADIV2S<-function(X,y,Z,trace=0,criteria="AICc",...){
	#---------------------------------------------------
	if(nrow(X)!=length(y) | nrow(X)!=nrow(Z)) stop("number of samples 
			of X, y and Z are not the same!")
	#---------------------------------------------------	
	n=nrow(X)
	p=ncol(X)
	q=ncol(Z)
	#---------------------------------------------------
	Gh=matrix(nrow=(q+1),ncol=p)
	for(j in 1:p){
		model2=grpreg(Z, X[,j], nlambda=q, alpha=1, intercept=FALSE, 
					group=1:q, penalty="grSCAD", gamma=3.7, family="gaussian",
					group.multiplier=rep(1,q))
		gamhat=coef(model2)
		gamhat=gamhat[,-1]
		Xhats=cbind(1,Z)%*%gamhat
		resids=X[,j]-Xhats
		RSSs=apply(resids,2,function(x){sum(x^2)})
		dfs=apply(gamhat,2,function(x){sum(x != 0)})
		if(criteria=="AICc")
			CRTs=2*dfs+n*log(RSSs/n)+(2*dfs^2+2*dfs)/(n-dfs-1)
		if(criteria=="AIC")
			CRTs=2*dfs+n*log(RSSs/n)
		if(criteria=="BIC")
			CRTs=log(n)*dfs+n*log(RSSs/n)
		if(criteria=="BICc")
			CRTs=log(q)*dfs+n*log(RSSs/n)
		path = which.min(CRTs)
		Gh[,j]=gamhat[,path]
	}
	Xhat=cbind(1,Z)%*%Gh
	model2=grpreg(Xhat, y, nlambda=p, alpha=1, intercept=FALSE, 
				group=1:p, penalty="grSCAD", gamma=3.7, family="gaussian",
				group.multiplier=rep(1,p))
	bethat=coef(model2)
	yhats=cbind(1,Xhat)%*%bethat
	resids=y-yhats
	RSSs=apply(resids,2,function(x){sum(x^2)})
	dfs=apply(bethat,2,function(x){sum(x != 0)})
		if(criteria=="AICc")
			CRTs=2*dfs+n*log(RSSs/n)+(2*dfs^2+2*dfs)/(n-dfs-1)
		if(criteria=="AIC")
			CRTs=2*dfs+n*log(RSSs/n)
		if(criteria=="BIC")
			CRTs=log(n)*dfs+n*log(RSSs/n)
		if(criteria=="BICc")
			CRTs=log(p)*dfs+n*log(RSSs/n)
	path = which.min(CRTs)
	bh=bethat[,path]
	sel = Gh!=0
	out=list(coefficients=bh,selected=sel,gamma=Gh)
	return(out)
}