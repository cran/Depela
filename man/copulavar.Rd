\name{copulavar}
\alias{copulavar}
\title{ Copula-VAR Estimation }
\description{
  Estimate the VAR system when the innorvation terms are generated from some copula models. This function can also be used to estimate the copula model when the observations are dependent acroos the time.
}
\usage{
copulavar(z, oder = 0, dens, init)
}
\arguments{
  \item{z}{ observations, a N*K matrix }
  \item{oder}{ the order of the VAR system. If order=0, then the program will select the order automatically based on the PACF criterion. }
  \item{dens}{ user-defined function returning the copula density given the value of dependence parameter. The definition of dens should follow the rule that the first argument is the value of the dependence parameter, and the second argument is the observation vector. }
  \item{init}{ initial guess of the value of the dependence parameter. }
}
\value{
  \item{m.est}{A list containing the marginal parameter}
  \item{mledp}{The ML estimate of the dependence parameter}
}
\examples{
##
	# generate random samples of a copula-var system
	x<-gendepela(500,2,list(c(0.2,0.1),c(0.2,0.11)),"clayton",8,rep(1,2))
	# define dens function
	dens<-function(theta,dat){
		u<-dat[1]
		v<-dat[2]
		return((1+theta)*u^(-1-theta)*v^(-1-theta)*(-1+u^(-theta)+v^(-theta))^(-2-1/theta))
	}
	out<-copulavar(x,2,dens,3)
}