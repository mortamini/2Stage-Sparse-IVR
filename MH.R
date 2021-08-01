AMH<-function(X,y,Z,init = list(model)){
	p = ncol(X)
	n = length(y)
	q = ncol(Z)
	bh0 = init$model$coefficients
	Gh0 = init$model$gamma
	dfb0 = sum(bh0 != 0)
	dfG0 = apply(Gh0,2,function(x){sum(x != 0)})
	nu0=sum(bh0^2)/dfb0
	omega0=sum(Gh0^2)/dfG0
	p0 = dfb0/p
	pi0=sum(dfG0)/p/q
	#
	Z1=cbind(1,Z)
	X1=cbind(1,X)
	Xhat0=Z1%*%Gh0
	yhat0=X1%*%bh0
	resy0=y-yhat0
	resX0=X-Xhat0
	RSSy0=sum(resy0^2)
	RSSX0=t(resX0)%*%resX0 
	RSSXy0 = as.vector(t(resy0)%*%resX0)
	Sig = rbind(c(RSSy0,RSSXy0) , cbind(RSSXy0,RSSX0))/n
	#
	mu.beta=rep(0,p)
	nu0=rep(nu0,p)
	mu.gamma=rep(0,p*q)
	omega0=rep(omega0,p*q)
	Iterations=100000
	#
	prop.s=rep(1,2*p+2*p*q)
	#
	#theta = matrix(nrow=Iterations,ncol=2*p+2*p*q)
	#
	loglike<-function(theta){
		beta = theta[1:p]
		gamma = theta[(p+1):(p+p*q)]
		Gam = matrix(gamma,q,p)
		B = cbind(c(1,beta),rbind(rep(0,p),diag(p)))
		Om = B%*%Sig%*%t(B)
		loglike = 0
		for(i in 1:n){
			yx = c(y[i],X[i,])
			zgb = c(t(Z[i,])%*%Gam%*%beta,t(Z[i,])%*%Gam)
			dif = yx - zgb
			pow = as.vector(t(dif)%*%tcinv(Om)%*%dif)
			loglike = loglike - 0.5*pow
		}
		loglike
	}
	#
	logprior<-function(theta){
		beta = theta[1:p]
		gamma = theta[(p+1):(p+p*q)]
		eta = theta[(p+p*q):(2*p+p*q)]
		tta = theta[(2*p+p*q):(2*p+2*p*q)]
		eta = exp(eta)/(1+exp(eta))
		tta = exp(tta)/(1+exp(tta))
		tta = (tta>0.5)
		eta = (eta>0.5)
		logprior = 0
		logprior = logprior + log(prod(eta*dnorm(beta,0,nu0)+
					(1-eta)*(beta==0)))
		logprior = logprior + log(prod(tta*dnorm(gamma,0,omega0)+
					(1-tta)*(gamma==0)))
		logprior = logprior + sum(dbinom(eta,1,p0,log=TRUE))
		logprior = logprior + sum(dbinom(tta,1,pi0,log=TRUE))
		logprior 
	}
	#
	acc.prob=0
	theta.hat=0
	current.theta=rep(0,2*p+2*p*q)
	for(t in 1:Iterations){
		progress(t,Iterations)
		if(t>10 & t<=1000){
			acp = acc.prob/t
			if(acp<0.2) prop.s = prop.s * 0.95
			if(acp>0.3) prop.s = prop.s * 1.05
		}
		prop.theta<-rnorm(2*p+2*p*q,current.theta,prop.s)
		loga<-loglike(prop.theta)-loglike(current.theta)+
			logprior(prop.theta)-logprior(current.theta)
		u=runif(1)
		lu=log(u)
		if(lu<loga){
			current.theta=prop.theta
			acc.prob=acc.prob+1
		}
		if(t>50000)
			theta.hat=theta.hat+current.theta/(Iterations-50000)
	}
	acp = acc.prob/Iterations
	beta.hat = theta.hat[1:p]
	gamma.hat = matrix(theta.hat[(p+1):(p+p*q)],q,p)
	eta.hat = theta.hat[(p+p*q):(2*p+p*q)]
	tta.hat = theta.hat[(2*p+p*q):(2*p+2*p*q)]
	eta.hat = exp(eta.hat)/(1+exp(eta.hat))
	tta.hat = exp(tta.hat)/(1+exp(tta.hat))
	tta.hat = (tta.hat>0.5)
	eta.hat = (eta.hat>0.5)
	selected = matrix(tta.hat,q,p)
	output = list(acp =acp , coefficients =beta.hat , gamma=gamma.hat,
		eta.hat=eta.hat, selected = selected)
	output
}
