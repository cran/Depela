copulavar<-function(z,oder=0,dens,init){
	# purpose: Estimate the VAR system when the innorvation terms are generated from some copula models. This function can also be used to estimate the copula model when the observations are dependent acroos the time.
	# arguments: 
	#	z	-- observations, a n*k matrix
	#	oder	-- the order of the VAR system. If order=0, then the program will select the order automatically based on the PACF criterion.
	#	dens	-- user-defined function returning the copula density given the value of dependence parameter. The definition of dens should follow the rule that the first argument is the value of the dependence parameter, and the second argument is the observation vector.
	#	init	-- initial guess of the value of the dependence parameter. a K*1 vector.
	# output:
	#	A list of MLE estimate reporting the estimates of the marginal and dependence parameters.
	
	# Before get-started.
	library(stats4) # Load stats4.
	K<-ncol(z)
	N<-nrow(z)
	ucdf<-matrix(nc=K,nr=N)
	m.est<-list() # list containing the marginal parameter
	
	# Define a inner function returning the order of AR given the pacf output.
	pacf2order<-function(pacf,N){
		P<-length(pacf)
		order<-0
		for(i in 1:P){
			if (pacf[i]>=(2/sqrt(N-i))){
				order<-order+1
			} else {
				break
			}
		}
		return(order)
	}

	# 1st step: Estimate the marginal parameter
	if(oder==0){
		for(k in 1:K){
			x<-z[,k]
			pacf.x<-as.vector(pacf(x)$acf) # deterimine the order for each univariate time series.
			od<-pacf2order(pacf.x,N)
			fit.x<-arima(x,order=c(od,0,0),include.mean=TRUE) # fit an AR
			m.est[[k]]<-fit.x #save the estimation results of fitting marginal time series
			res.x<-fit.x$residuals # pick out the residuals

			cdfk<-ecdf(res.x) # transform to their cdf level
			ucdf[,k]<-cdfk(res.x)
		}
	}else{
		for(k in 1:K){
			x<-z[,k]
			fit.x<-arima(x,order=c(oder,0,0),include.mean=TRUE) # fit an AR
			m.est[[k]]<-fit.x #save the estimation results of fitting marginal time series
			res.x<-fit.x$residuals # pick out the residuals

			cdfk<-ecdf(res.x) # transform to their cdf level
			ucdf[,k]<-cdfk(res.x)
		}
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

	# report estimation result
	out.est<-list(m.est=m.est,mledp=mledp)
	return(out.est)
}