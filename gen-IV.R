gen.IV<-function(beta,Gam,Z,s02=1,tau02=1,int1=0,int2=0){
	#
	y=c()
	X=matrix(nrow=n,ncol=p)
	for(i in 1:n){
		for(j in 1:p){
			X[i,j]=rnorm(1,int1+as.vector((Z[i,])%*%Gam[,j]),tau02)		
		}
	y[i]=rnorm(1,int2+t(X[i,])%*%beta,s02)
	}
	out=list(X=X,y=y)
}