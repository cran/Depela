\name{spestimator}
\alias{spestimator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two-Stage Semiparametric Estimator of Copula Model }
\description{
  \code{spestimator} estimates the dependence parameter of a copula model given observations by implementing the semi-parametric estimator. 
}
\usage{
spestimator(u, dens, init)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{ observations, a n*k matrix }
  \item{dens}{ user-defined function returning the copula density given the value of dependence parameter. The definition of dens should follow the rule that the first argument is the value of the dependence parameter, and the second argument is the observation vector. }
  \item{init}{ initial guess of the value of the dependence parameter. a k*1 vector. }
}
\details{
  \code{spestimator} first estimate the marginal CDF for each random variable using nonparametric method, then estimate the dependence parameter by implementing MLE.
}
\value{
  A ML estimate of the dependence parameter which is an object of class mle-class.
}

\examples{
##	
	# Example: Two-stage estimation of a two-dimensional Clayton copula model.
	## Define the density function for Clayton copula.
	library(copula)
	dens<-function(theta,dat){
		u<-dat[1]
		v<-dat[2]
		return((1+theta)*u^(-1-theta)*v^(-1-theta)*(-1+u^(-theta)+v^(-theta))^(-2-1/theta))
	}
	x <- mvdc(claytonCopula(8), c("norm", "exp"),list(list(mean = 0, sd =2), list(rate = 2)))
	x.samp <- rmvdc(x, 500)
	tsmle<-spestimator(x.samp,dens,init=2)
}