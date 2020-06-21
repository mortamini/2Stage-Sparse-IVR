cv.model <- function(model.name,fold=3,...){
	model <- match.fun(model.name)
	cvc=c()
	nts = trunc(n/fold)
	for(f in 1:fold){
		test= nts*(f-1)+1:nts
		train.X = X[-test,]
		train.Z = Z[-test,]
		train.y = y[-test]
		#
		test.Z = Z[test,]
		test.X = X[test,]
		test.y = y[test]
		fit<-model(train.X,train.y,train.Z,...)
		bh = fit$coefficients
		Gh = fit$gamma
		p=ncol(test.X)
		q=ncol(test.Z)
		Z1=test.Z
		if(nrow(Gh)==(q+1)) Z1=cbind(1,test.Z) 
		Xhat=Z1%*%Gh
		if(length(bh)==(p+1)) Xhat=cbind(1,Xhat) 
		yhat=Xhat%*%bh
		cvc[f] = sqrt(sum(yhat - test.y)^2/n)
	}
	cv = mean(cvc)
	return(cv)
}