\name{gendepela}
\alias{gendepela}
\title{ Generate Copula VAR Time Series }
\description{
  \code{gendepela} generate multivariate time series allowing significant richer dependence structure among innovation terms.
}
\usage{
gendepela(N, d, phi, copfamily, theta, bstd)
}
\arguments{
  \item{N}{ number of observations }
  \item{d}{ the order of the VAR system }
  \item{phi}{ a list containing the value of marginal parameter. The length of each element of this list should be equal }
  \item{copfamily}{ the family name of the copula. "clayton", "frank", and "gumbel" are supported }
  \item{theta}{ the value of dependence parameter }
  \item{bstd}{ a vector containing the values of sd of for each innovation term }
}
\details{
  The marginal distribution of each disturbance term is normal distribution with mean 0, and sd assgined by bstd. However, their dependence is determined by the assgined copula parameterized by theta.
}
\value{
  A data set generated from user-defined depela model.
}
\examples{
##
	# Generate a 50*3 multivariate time series whose innovation terms are contructed from a Gumbel copula.
	 dat<-gendepela(50,2,list(c(0.2,0.1),c(0.2,0.11),c(0.2,0.1)),"gumbel",4,rep(1,3))
}