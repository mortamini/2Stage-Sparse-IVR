EPIV2S <- function(X, y, Z, trace=2, eps1=1e-5, eps2=1e-5, 
			post=FALSE, intercept = TRUE,maxit1=20,maxit2=30,
			init=list(p0=0.5,pi0=0.5,v02=1,w02=1,s02=1,tau02=1)){
	# init = list(p0,pi0,v02,w02,s02,tau02)
				# cross-validation or evaluation on all parameters
	# init = list(model,p0,pi0)
				# Strategy I - cross-validation or evaluation on p0,pi0
	# init = list(model,v02,w02)
				# Strategy II - cross-validation or evaluation on v02,w02

	#---------------------------------------------------
	if(nrow(X)!=length(y) | nrow(X)!=nrow(Z)) stop("number of samples 
			of X, y and Z are not the same!")
	#---------------------------------------------------	
	n=nrow(X)
	p=ncol(X)
	q=ncol(Z)
	#---------------------------------------------------
	p2=p2_old=v2=v2_old=m2=m2_old=v1=v1_old=m1=m1_old=c()
	pi2=pi2_old=w2=w2_old=mu2=mu2_old=w1=w1_old=mu1=mu1_old=c()
	#---------------------------------------------------
	w1=rep(1e+100,p*q)
	w2=rep(1e+100,p*q)
	#---------------------------------------------------
	if(is.null(init$model)){
		p0=init$p0
		pi0=init$pi0
		v02=init$v02
		w02=init$w02
		s02=init$s02
		tau02=init$tau02
	}
	if(is.null(init$v02) & is.null(init$w02)){
		p0=init$p0
		pi0=init$pi0
		v02=1
		w02=1
		bh0 = init$model$coefficients
		Gh0 = init$model$gamma
		Z1=cbind(1,Z)
		X1=cbind(1,X)
		Xhat0=Z1%*%Gh0
		yhat0=X1%*%bh0
		resy0=y-yhat0
		resX0=X-Xhat0
		RSSy0=sum(resy0^2)
		RSSX0=sum(resX0^2)
		dfye=ifelse(dfb0>=n,trunc(n/2),n-dfb0)
		dfXe=ifelse(sum(dfG0)>=n*p,p*trunc(n/2),p*n-sum(dfG0))
		s02=min(1,max(0.01,RSSy0/dfye/p0))
		tau02=min(1,max(0.01,RSSX0/n))
	}
	if(is.null(init$p0) & is.null(init$pi0)){
		v02=init$v02
		w02=init$w02
		bh0 = init$model$coefficients
		Gh0 = init$model$gamma
		dfb0 = sum(bh0 != 0)
		dfG0 = apply(Gh0,2,function(x){sum(x != 0)})
		p0 = dfb0/p
		pi0=sum(dfG0)/p/q
		Z1=cbind(1,Z)
		X1=cbind(1,X)
		Xhat0=Z1%*%Gh0
		yhat0=X1%*%bh0
		resy0=y-yhat0
		resX0=X-Xhat0
		RSSy0=sum(resy0^2)
		RSSX0=sum(resX0^2)
		dfye=ifelse(dfb0>=n,trunc(n/2),n-dfb0)
		dfXe=ifelse(sum(dfG0)>=n*p,p*trunc(n/2),p*n-sum(dfG0))
		s02=min(1,max(0.01,RSSy0/dfye/p0))
		tau02=min(1,max(0.01,RSSX0/n))
	}
	#---------------------------------------------------
	p3=log(p0/(1-p0))
	pi3=log(pi0/(1-pi0))
	#---------------------------------------------------
	if(intercept){
		Zm=apply(Z,2,mean)
		Xm=apply(X,2,mean)
		Znc=Z
		Xnc=X
		Zc=t(t(Z)-Zm)
		Xc=t(t(X)-Xm)
		Z=Zc
		X=Xc
	}
	p=ncol(X)
	q=ncol(Z)
	mu1=rep(0,p*q)
	mu2=rep(0,p*q)
	pi2=rep(0,p*q)
	#---------------------------------------------------
	s=function(H){
		ms=1/(1+exp(-H))
		return(ms)
	}
	#---------------------------------------------------
	mu1norm=mu2norm=w1norm=w2norm=1
	iter=1
	e=1
	red=0.99
	#---------------------------------------------------
	if(trace >= 1) cat("Running 2-stage EP algorithm ...\n")
	if(trace >= 1) cat("Stage I ...\n")
	err1 = (mu1norm+mu2norm+w1norm+w2norm)
	if(trace == 2) cat("EP iteration      error \n")
	while(err1 >= eps1 & iter <=maxit1){
		err1=(mu1norm+mu2norm+w1norm+w2norm)
		if(trace == 2) cat(paste(iter,"                   ",err1,"\n"))
		#---------------------------------------------------
		iter=iter+1
		if(iter>=10) red = red^2
		#---------------------------------------------------
		pi2_old= pi2
		w1_old=w1
		w2_old=w2
		mu1_old=mu1
		mu2_old=mu2
		#---------------------------------------------------
		pi2 = 1/2 * log(w1) - 1/2 * log(w1+w02) + 1/2 * mu1^2 * (1/w1-1/(w1+w02)) 
		pi2_damp = e * pi2 + (1-e) * pi2_old
		#---------------------------------------------------
		c = s(pi2+pi3) * mu1/(w1+w02) + s(-(pi2+pi3))* mu1/w1
		d = s(pi2+pi3) * (mu1^2-w1-w02)/(w1+w02)^2 + s(-(pi2+pi3)) * (mu1^2/w1^2 - 1/w1)
		#---------------------------------------------------
		w2 = 1/(c^2 - d) - w1
		for(j in 1:(p*q)){
			if(w2[j]<=0 || w2[j]>100 ) w2[j]=100
		}
		w2_damp=1/(e/w2 + (1-e)/w2_old)
		#---------------------------------------------------
		mu2 = mu1 - c* (w2 + w1)
		mu2_damp = w2_damp * (e * mu2/w2 + (1-e) * mu2_old/w2_old)
		#---------------------------------------------------
		mu2=mu2_damp
		w2=w2_damp
		pi2=pi2_damp
		#---------------------------------------------------
		Q1=t(Z)%*%X
		wnew=c()
		munew=c()
		for(i in 1:p){
			progress(x=i,max=p)
			mu2i=mu2[q*(i-1)+(1:q)]
			w2i=diag(w2[q*(i-1)+(1:q)])
			W2z=Z%*%w2i
			Wi=w2i-t(W2z)%*%solve(tau02*diag(n)+W2z%*%t(Z))%*%W2z
			munew=c(munew,Wi%*%(mu2i/diag(w2i)+tau02^(-1)*Q1[,i]))
			#
			wnew=c(wnew,diag(Wi))
		}
		#
		#---------------------------------------------------
		w1=1/(1/wnew-1/w2)
		w1_damp=1/(e/w1 + (1-e)/w1_old)
		#---------------------------------------------------
		mu1=(munew/wnew-mu2/w2)*w1
		mu1_damp = w1_damp * (e * mu1/w1 + (1-e) * mu1_old/w1_old)
		#---------------------------------------------------
		mu1=mu1_damp
		w1=w1_damp
		#---------------------------------------------------
		mu1norm=as.vector(t(mu1-mu1_old)%*%(mu1-mu1_old))
		mu2norm=as.vector(t(mu2-mu2_old)%*%(mu2-mu2_old))
		w1norm=as.vector(t(w1-w1_old)%*%(w1-w1_old))
		w2norm=as.vector(t(w2-w2_old)%*%(w2-w2_old))
		#---------------------------------------------------
		#
		e=e*red
		#---------------------------------------------------
	}
	w1=as.vector(w1)
	w2=as.vector(w2)
	Gammahat=matrix(1/(1/w1+1/w2)*(mu1/w1+mu2/w2),q,p)
	Xhat=Z%*%Gammahat
	if(intercept){
		Gamma0=as.vector(Xm-Zm%*%Gammahat)
		Xhat=t(Gamma0+t(Z%*%Gammahat))
		Gammahat=rbind(Gamma0,Gammahat)
		Xhatm=apply(Xhat,2,mean)
		Xhatnc=Xhat
		Xhatc=t(t(Xhat)-Xhatm)
		Xhat=Xhatc
		ym=mean(y)
		yc=y-ym
		ync=y
		y=yc
	}
	#---------------------------------------------------
	v1=rep(1e+100,p)
	v2=rep(1e+100,p)
	m1=rep(0,p)
	m2=rep(0,p)
	p2=rep(0,p)
	e=1
	red=0.99
	#---------------------------------------------------
	m1norm=m2norm=v1norm=v2norm=1
	iter=0
	#---------------------------------------------------
	if(trace >= 1) cat("Stage II ...\n")
	if(trace == 2) cat("EP iteration      error \n")
	err2=(m1norm+m2norm+v1norm+v2norm)
	while(err2 >= eps2 & iter<=maxit2){
		#---------------------------------------------------
		err2=(m1norm+m2norm+v1norm+v2norm)
		if(trace == 2) cat(paste(iter,"      ",err2,"\n"))
		#---------------------------------------------------
		iter=iter+1
		if(iter>=10) red = red^2
		#---------------------------------------------------
		p2_old= p2
		v1_old=v1
		v2_old=v2
		m1_old=m1
		m2_old=m2
		#---------------------------------------------------
		p2 = 1/2 * log(v1) - 1/2 * log(v1+v02) + 1/2 * m1^2 * (1/v1-1/(v1+v02)) 
		p2_damp = e * p2 + (1-e) * p2_old
		#---------------------------------------------------
		a = s(p2+p3) * m1/(v1+v02) + s(-(p2+p3))* m1/v1
		b = s(p2+p3) * (m1^2-v1-v02)/(v1+v02)^2 + s(-(p2+p3)) * (m1^2/v1^2 - 1/v1)
		#---------------------------------------------------
		v2 = 1/(a^2 - b) - v1
		for(j in 1:p){
			if(v2[j]<=0 || v2[j]>100 ) v2[j]=100
		}
		v2_damp=1/(e/v2 + (1-e)/v2_old)
		#---------------------------------------------------
		m2 = m1 - a* (v2 + v1)
		m2_damp = v2_damp * (e * m2/v2 + (1-e) * m2_old/v2_old)
		#---------------------------------------------------
		m2=m2_damp
		v2=v2_damp
		p2=p2_damp
		#---------------------------------------------------
		V2=diag(as.vector(v2))
		V2x=V2%*%t(Xhat)
		V=V2-V2x%*%ginv(s02*diag(n)+Xhat%*%V2x)%*%t(V2x)
		mnew=V%*%(m2/v2+s02^(-1)*t(Xhat)%*%y)
		vnew=diag(V)
		#---------------------------------------------------
		#
		v1=1/(1/vnew-1/v2)
		v1_damp=1/(e/v1 + (1-e)/v1_old)
		#---------------------------------------------------
		m1=(mnew/vnew-m2/v2)*v1
		m1_damp = v1_damp * (e * m1/v1 + (1-e) * m1_old/v1_old)
		#---------------------------------------------------
		m1=m1_damp
		v1=v1_damp
		#---------------------------------------------------
		m1norm=as.vector(t(m1-m1_old)%*%(m1-m1_old))
		m2norm=as.vector(t(m2-m2_old)%*%(m2-m2_old))
		v1norm=as.vector(t(v1-v1_old)%*%(v1-v1_old))
		v2norm=as.vector(t(v2-v2_old)%*%(v2-v2_old))
		#---------------------------------------------------
		e=e*red
		#---------------------------------------------------
	}
	#---------------------------------------------------
	v1=as.vector(v1)
	v2=as.vector(v2)
	#---------------------------------------------------
	betahat=1/(1/v1+1/v2)*(m1/v1+m2/v2)
	betasel=s(-p2-p3)<=quantile(s(-p2-p3),p0)
	gammasel=matrix(s(-pi2-pi3)<=quantile(s(-pi2-pi3),pi0),q,p)
	if(intercept & !post){
		beta0=ym-Xhatm%*%betahat
		betahat=c(beta0,betahat)
		betasel=c(TRUE,betasel)
		gammasel=rbind(TRUE,gammasel)
	}
	var1=1/(1/v1+1/v2)
	var2=1/(1/w1+1/w2)
	if(post){
		if(trace>=1) cat("Post estimation ...\n")
		Gammahat=SigGam=matrix(0,nrow=q,ncol=p)
		for(j in 1:p){
			ZSj=Z[,gammasel[,j]]
			if(sum(gammasel[,j])>0 & sum(gammasel[,j])<n){
				SigGamj=tau02*ginv(t(ZSj)%*%ZSj)
				Gammahat[gammasel[,j],j]=SigGamj%*%t(ZSj)%*%X[,j]/tau02
				SigGam[gammasel[,j],j]=diag(SigGamj)
			}else if(sum(gammasel[,j])>=n){
				SigGamj=ginv(t(ZSj)%*%ZSj/tau02+diag(sum(gammasel[,j]))/w02)
				Gammahat[gammasel[,j],j]=SigGamj%*%t(ZSj)%*%X[,j]/tau02
				SigGam[gammasel[,j],j]=diag(SigGamj)
			}
		}
		if(intercept){
			Gamma0=as.vector(Xm-Zm%*%Gammahat)
			Xhat=t(Gamma0+t(Z%*%Gammahat))
			Xhatm=apply(Xhat,2,mean)
			Gammahat=rbind(Gamma0,Gammahat)
			Xhatc=t(t(Xhat)-Xhatm)
			Xhat=Xhatc
		}else{
			Xhat=Z%*%Gammahat
		}
		if(sum(betasel)<n){
			XhatS=Xhat[,betasel]
			SigBet=s02*ginv(t(XhatS)%*%XhatS)
			betahat[betasel]=SigBet%*%t(XhatS)%*%y/s02
		}else{
			XhatS=Xhat[,betasel]
			SigBet=ginv(t(XhatS)%*%XhatS/s02+diag(sum(betasel))/v02)
			betahat[betasel]=SigBet%*%t(XhatS)%*%y/s02
		}
		betahat[!betasel]=0
		if(intercept){
			beta0=ym-Xhatm%*%betahat
			betahat=c(beta0,betahat)
		}
		se1=rep(0,p)
		se1[betasel]=sqrt(diag(SigBet))
		se2=matrix(sqrt(as.vector(SigGam)),q,p)
	}else{
		betahat[!betasel]=0
		Gammahat[!gammasel]=0
		se1=sqrt(var1)
		se2=matrix(sqrt(var2),q,p)
	}
	#
	if(intercept & post) 	gammasel=rbind(TRUE,gammasel)
	out=list(coefficients=betahat,selected=gammasel,gamma=Gammahat,
		se1=se1,se2=se2,err1=err1,err2=err2)
	return(out)
}
