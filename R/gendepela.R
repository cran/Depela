gendepela<-function(N,d,phi,copfamily,theta,bstd){
	# purpose: Generate copula VAR time series.
	# arguments: 
	#	N	-- number of observations
	#	d	-- the order of the VAR system
	#	phi	-- a list containing the value of marginal parameter. The length of each element of this list should be equal
	#	copfamily	-- the family name of the copula. "clayton", "frank", and "gumbel" are supported.
	#	theta	-- the value of dependence parameter
	#	bstd	-- a vector containing the values of sd of for each innovation term.
	# details:
	#	The marginal distribution of each disturbance term is normal distribution with mean 0, and sd assgined by bstd. However, their dependence is determined by the assgined copula parameterized by theta.
	# output:
	#	A data set.

	# preparement
	K<-length(phi) # the rank of the VAR
	
	# copula settings
	library(copula)
	clay<-archmCopula(family=copfamily,dim=K,param=theta) # Generate a user-defined copula

	# Data generation.
	# define the normal margins
	normm<-c()
	pm<-list()
	for(k in 1:K){
		normm[k]<-"norm"
		pm[[k]]<-list(mean=0,sd=bstd[k])
	}

	mvd<-mvdc(copula=clay,margins=normm,paramMargins=pm)
	dat<-rmvdc(mvd,(N+d)) # "dat" is the error term series [epsilon1(t),epsilon2(t)] which copuled by a Clayton copula with N(0,1) margins.
	y<-matrix(NA,nr=N+d,nc=K)

	for(k in 1:K){
		y[1:d,k]<-rep(0,d) # initial value
	}

	for (i in (d+1):(N+d)){
		for(k in 1:K){
			y[i,k]<-dat[i,k]+sum(phi[[k]]*y[(i-d[1]):(i-1),k]) # modification wrt d
		}
	}
	y<-y[(d+1):(N+d),]
	return(y)
}