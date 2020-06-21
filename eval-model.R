eval.model<-function(model){
	bh=model$coefficients
	sel=model$selected
	Gh=model$gamma
	p=ncol(X)
	q=ncol(Z)
	Z1=Z
	if(nrow(Gh)==(q+1)) Z1=cbind(1,Z) 
	Xhat=Z1%*%Gh
	if(length(bh)==(p+1)) Xhat=cbind(1,Xhat) 
	yhat=Xhat%*%bh
	RMSPE=sqrt(sum((yhat-y)^2)/n)
	if(nrow(Gh)==(q+1)) sel=sel[-1,]
	if(length(bh)==(p+1)) bh=bh[-1]
	sum11=sum(beta==0 & bh!=0)
	sum21=sum(beta==0)
	FPRy=(sum11/sum21)
	#
	sum12=sum(Gam==0 & sel)
	sum22=sum(Gam==0)
	FPRX=(sum12/sum22)
	#
	sum13=sum(beta!=0 & bh==0)
	sum23=sum(beta!=0)
	FNRy=(sum13/sum23)
	#
	sum14=sum(Gam!=0 & !sel)
	sum24=sum(Gam!=0)
	FNRX=(sum14/sum24)
	out=c(FPRy,FPRX,FNRy,FNRX,RMSPE)
	return(out)
}