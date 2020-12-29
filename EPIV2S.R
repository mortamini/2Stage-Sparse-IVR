EPIV2S <- function(X, y, Z, trace=2, eps=c(1e-5, 1e-5), 
			post=FALSE, intercept = c(TRUE,TRUE),maxit=c(20,30),
			init=list(p0=0.5,pi0=0.5,v02=1,w02=1,s02=1,tau02=1)){
	# init = list(p0,pi0,v02,w02,s02,tau02)
				# cross-validation or evaluation on all parameters
	# init = list(model,p0,pi0)
				# cross-validation or evaluation on p0,pi0
	# init = list(model)
				# initial all from model 

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
	scalex=sd(X)
	scalez=sd(Z)
	scaley=sd(y)
	#---------------------------------------------------
	if(is.null(init$model)){
		p0=init$p0
		pi0=init$pi0
		v02=init$v02
		w02=init$w02
		s02=init$s02
		tau02=init$tau02
	}
	if(!is.null(init$model) & !is.null(init$p0) & !is.null(init$pi0)){
		p0=init$p0
		pi0=init$pi0
		bh0 = init$model$coefficients
		Gh0 = init$model$gamma
		dfb0 = sum(bh0 != 0)
		dfG0 = apply(Gh0,2,function(x){sum(x != 0)})
		v02=sum(bh0^2)/dfb0
		w02=sum(Gh0^2)/dfG0
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
		s02=RSSy0/dfye/scaley^2
		tau02=RSSX0/n/scalex^2
#		s02=v02*RSSy0/dfye/scaley^2/sum(bh0^2)*dfb0
#		tau02=w02*RSSX0/dfXe/scalex^2/sum(Gh0%*%t(Gh0))*sum(dfG0)
#		cat("tau02= \n")
#		print(tau02)
#		cat("s02= \n")
#		print(s02)
		s02=min(scaley^2,max(0.01,s02))
		tau02=min(scalex^2,max(0.01,tau02))
	}
	if((!is.null(init$model) & is.null(init$p0) & is.null(init$pi0))){
		bh0 = init$model$coefficients
		Gh0 = init$model$gamma
		dfb0 = sum(bh0 != 0)
		dfG0 = apply(Gh0,2,function(x){sum(x != 0)})
		v02=sum(bh0^2)/dfb0
		w02=sum(Gh0^2)/dfG0
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
		s02=RSSy0/dfye/scaley^2
		tau02=RSSX0/n/scalex^2
#		s02=v02*RSSy0/dfye/scaley^2/sum(bh0^2)*dfb0
#		tau02=w02*RSSX0/dfXe/scalex^2/sum(Gh0%*%t(Gh0))*sum(dfG0)
#		cat("tau02= \n")
#		print(tau02)
#		cat("s02= \n")
#		print(s02)
		s02=min(scaley^2/p0,max(0.01,s02))
		tau02=min(scalex^2/pi0,max(0.01,tau02))
	}
	#---------------------------------------------------
	p3=log(p0/(1-p0))
	pi3=log(pi0/(1-pi0))
	#---------------------------------------------------
		X=X/scalex
		Z=Z/scalez
		y=y/scaley
	#---------------------------------------------------
	if(intercept[1]){
		Xm=apply(X,2,mean)
		Xnc=X
		Xc=t(t(X)-Xm)
		X=Xc
	}
	if(intercept[2]){
		Zm=apply(Z,2,mean)
		Znc=Z
		Zc=t(t(Z)-Zm)
		Z=Zc
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
	al0=3*pi0
	be0=3*(1-pi0)
	aa0=3/2
	bb0=3/2/tau02
	a0p=3/2
	b0p=3/2/w02
	#---------------------------------------------------
	if(trace >= 1) cat("Running 2-stage EP algorithm ...\n")
	if(trace >= 1) cat("Stage I ...\n")
	err1 = (mu1norm+mu2norm+w1norm+w2norm)
	if(trace == 2) cat("EP iteration      error \n")
	while(err1 >= eps[1] & iter <=maxit[1]){
		err1=(mu1norm+mu2norm+w1norm+w2norm)
		if(trace == 2) cat(paste(iter,"                   ",err1,"\n"))
		#---------------------------------------------------
		iter=iter+1
		if(iter>=10) red = red^2
		#---------------------------------------------------
		if(iter >= 3){
			sz = sum(s(pi2+pi3))
#			pi0 = (al0+sz)/(al0+be0+p*q)
#			pi0 = max(0.05,pi0)
			Gh0=matrix(1/(1/w1+1/w2)*(mu1/w1+mu2/w2),q,p)
			gsel0=matrix(s(-pi2-pi3)<=quantile(s(-pi2-pi3),pi0),q,p)
			if(intercept[2]){
				G0=as.vector(Xm-Zm%*%Gh0)
				Xhat0=t(G0+t(Z%*%Gh0))
				Gh0=rbind(G0,Gh0)
				gsel0=rbind(TRUE,gsel0)
			}else{
				Xhat0=Z%*%Gh0
			}
			resX0=X-Xhat0
			RSSX0=sum(resX0^2)
#			tau02 = (bb0+RSSX0/2)/(aa0+n/2-1)
			Gz=sum(Gh0^2*gsel0)
			w02 = (b0p+Gz)/(a0p+sz/2-1)			
		}
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
#		cat("w2 before trim= \n")
#		print(summary(w2))
		for(j in 1:(p*q)){
			if(w2[j]<=0 || w2[j]>100 ) w2[j]=100
		}
#		cat("w2 after trim= \n")
#		print(summary(w2))
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
			if(trace>=2) progress(x=i,max=p)
			mu2i=mu2[q*(i-1)+(1:q)]
			w2i=diag(w2[q*(i-1)+(1:q)])
			W2z=Z%*%w2i
			Wi=w2i-t(W2z)%*%tcinv(tau02*diag(n)+W2z%*%t(Z))%*%W2z
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
#		cat("w1 before trim= \n")
#		print(summary(w1))
		for(j in 1:(p*q)){
			if(w1[j]<=0) w1[j]=1e-100
		}
#		cat("w1 after trim= \n")
#		print(summary(w1))
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
	cat("pio		wo2		tau02	\n")
	cat(pi0,"	",w02,"		",tau02,"	\n")
	w1=as.vector(w1)
	w2=as.vector(w2)
	Gammahat=matrix(1/(1/w1+1/w2)*(mu1/w1+mu2/w2),q,p)
	Xhat=Z%*%Gammahat
	if(intercept[2]){
		Gamma0=as.vector(Xm-Zm%*%Gammahat)
		Xhat=t(Gamma0+t(Z%*%Gammahat))
		Gammahat=rbind(Gamma0,Gammahat)
	}
	that = sum(X - Xhat)^2/n
	if(intercept[1]){
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
	al0=3*p0
	be0=3*(1-p0)
	aa0=3/2
	bb0=3/2/s02
	a0p=3/2
	b0p=3/2/v02
	#---------------------------------------------------
	m1norm=m2norm=v1norm=v2norm=1
	iter=0
	#---------------------------------------------------
	if(trace >= 1) cat("Stage II ...\n")
	if(trace == 2) cat("EP iteration      error \n")
	err2=(m1norm+m2norm+v1norm+v2norm)
	while(err2 >= eps[2] & iter<=maxit[2]){
		#---------------------------------------------------
		err2=(m1norm+m2norm+v1norm+v2norm)
		if(trace == 2) cat(paste(iter,"      ",err2,"\n"))
		#---------------------------------------------------
		iter=iter+1
		if(iter>=10) red = red^2
		#---------------------------------------------------
		if(iter >= 3){
			sz = sum(s(p2+p3))
#			p0 = (al0+sz)/(al0+be0+p)
#			p0 = max(0.05,p0)
			bet0=1/(1/v1+1/v2)*(m1/v1+m2/v2)
			bsel0=s(-p2-p3)<=quantile(s(-p2-p3),p0)
			if(intercept[1]){
				b0=as.vector(ym-Xhatm%*%bet0)
				yhat0=t(b0+t(Xhat%*%bet0))
				bet0=rbind(b0,bet0)
				bsel0=rbind(TRUE,bsel0)
			}else{
				yhat0=Xhat%*%bet0
			}
			resy0=y-yhat0
			RSSy0=sum(resy0^2)
#			s02 = (bb0+RSSy0/2)/(aa0+n/2-1)
			bz=sum(bet0^2*bsel0)
			v02 = (b0p+bz)/(a0p+sz/2-1)			
		}
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
		V=V2-V2x%*%tcinv(s02*diag(n)+Xhat%*%V2x)%*%t(V2x)
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
	cat("p0		v02		s02	\n")
	cat(p0,"	",v02,"		",s02,"	\n")
	#---------------------------------------------------
	betahat=1/(1/v1+1/v2)*(m1/v1+m2/v2)
	betasel=s(-p2-p3)<=quantile(s(-p2-p3),p0)
	if(intercept[1]){
		b0=as.vector(ym-Xhatm%*%betahat)
		yhat0=t(b0+t(Xhat%*%betahat))
		bet0=rbind(b0,betahat)
		bsel0=rbind(TRUE,betasel)
	}else{
		yhat0=Xhat%*%betasel
	}
	resy0=y-yhat0
	shat=sum(resy0^2)/n
	sb2 = sum(betahat^2)
	gammasel=matrix(s(-pi2-pi3)<=quantile(s(-pi2-pi3),pi0),q,p)
	if(intercept[1] & !post){
		beta0=ym-Xhatm%*%betahat
		betahat=c(beta0,betahat)
		betasel=c(TRUE,betasel)
	}
	if(intercept[2] & !post){
		gammasel=rbind(TRUE,gammasel)
	}
	var1=1/(1/v1+1/v2)
	var2=1/(1/w1+1/w2)
	if(post){
		if(trace>=1) cat("Post estimation ...\n")
		if(intercept[2]){Gh0=Gammahat[-1,]}else{Gh0=Gammahat}
		Gammahat=SigGam=matrix(0,nrow=q,ncol=p)
		for(j in 1:p){
			if(trace>=2) progress(x=j,max=p)
			ZSj=Z[,gammasel[,j]]
			if(sum(gammasel[,j])>0){
#			if(sum(gammasel[,j])>0 & sum(gammasel[,j])<n){
#				SigGamj=tcinv(t(ZSj)%*%ZSj)
#				Gammahat[gammasel[,j],j]=SigGamj%*%t(ZSj)%*%X[,j]
#				SigGam[gammasel[,j],j]=diag(SigGamj)
#			}else if(sum(gammasel[,j])>=n){
				sgam2 =sum(Gh0[,j]^2)
				lambda1= p * q * pi0* that / sgam2
				lambda1 = min(0.01,lambda1)
				SigGamj=tcinv(t(ZSj)%*%ZSj+lambda1*diag(sum(gammasel[,j])))
				Gammahat[gammasel[,j],j]=SigGamj%*%t(ZSj)%*%X[,j]
				SigGam[gammasel[,j],j]=diag(SigGamj)
#			}
			}
		}
		if(intercept[2]){
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
			SigBet=tcinv(t(XhatS)%*%XhatS)
			betahat[betasel]=SigBet%*%t(XhatS)%*%y
		}else{
			XhatS=Xhat[,betasel]
			lambda2 = p * p0 * shat /sb2
			lambda2 = min(0.01,lambda2)
			SigBet=tcinv(t(XhatS)%*%XhatS+lambda2*diag(sum(betasel)))
			betahat[betasel]=SigBet%*%t(XhatS)%*%y
		}
		betahat[!betasel]=0
		if(intercept[1]){
			beta0=ym-Xhatm%*%betahat
			betahat=c(beta0,betahat)
		}
		se1=rep(0,p)
		se1[betasel]=sqrt(shat*diag(SigBet))
		se2=matrix(sqrt(that*as.vector(SigGam)),q,p)
	}else{
		betahat[!betasel]=0
		Gammahat[!gammasel]=0
		se1=sqrt(var1)
		se2=matrix(sqrt(var2),q,p)
	}
	#
	if(intercept[2] & post) 	gammasel=rbind(TRUE,gammasel)
	if(intercept[1]){
		betahat=c(betahat[1]*scaley,betahat[-1]/scalex*scaley)
	}else{
		betahat=betahat/scalex*scaley
	}
	if(intercept[2]){
		Gammahat=cbind(Gammahat[,1]*scalex,Gammahat[,-1]/scalez*scalex)	
	}else{
		Gammahat=Gammahat/scalez*scalex
	}
	se1=se1/scalex*scaley
	se2=se2/scalez*scalex
	out=list(coefficients=betahat,selected=gammasel,gamma=Gammahat,
		se1=se1,se2=se2,err1=err1,err2=err2)
	return(out)
}
