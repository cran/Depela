spestimator<-function(u,dens,init){
	# purpose: Implement the semi-parametric estimator for estimating copula model given observations u.
	# arguments: 
	#	u	-- observations, a n*k matrix
	#	dens	-- user-defined function returning the copula density given the value of dependence parameter. The definition of dens should follow the rule that the first argument is the value of the dependence parameter, and the second argument is the observation vector.
	#	init	-- initial guess of the value of the dependence parameter. a K*1 vector.
	# output:
	# 	A ML estimate of the dependence parameter which is an object of class mle-class.
	
	# Before get-started.
	library(stats4) # Load stats4.
	K<-ncol(u)
	N<-nrow(u)
	ucdf<-matrix(nc=K,nr=N)
	
	# 1st step: calculate the empirical CDF for each random variable
	for(k in 1:K){
		cdfk<-ecdf(u[,k])
		ucdf[,k]<-cdfk(u[,k])
	}
	
	# Defien negative log-likelihood function
	like<-function(depv){
		lk<-0
		for(i in 1:N){
			lk<--log(dens(depv,ucdf[i,]))+lk
		}
		return(lk)#
	}
	# 2nd step: Implement the MLE to estimate the dependence parameter.
	mledp<-mle(like,start=list(depv=init))

	# report the ML estimate of the dependence parameter.
	return(mledp)
}