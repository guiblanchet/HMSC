\name{communitySimul}
\alias{communitySimul}

\title{
Simulate a species community
}
\description{
Simulates a species community assuming a linear model with the proposed descriptors, the species traits and one or more random effects.
}
\usage{
communitySimul(X,Tr=NULL,Random=NULL,nsp=NULL,family="probit",
			   paramX=NULL,paramTr=NULL,meansparamX=NULL,varX=NULL,
			   latent=NULL,paramLatent=NULL,spCor=NULL)
}
\arguments{
  \item{X}{
Matrix of descriptors defining the species. Rows are sites and the columns are descriptors.
}
  \item{Tr}{
Matrix defining species traits. Rows are traits and columns are species.
}
  \item{Random}{
A factor or a data.frame that includes many factors defining a random effect on the sites.
}
  \item{nsp}{
Numeric. Number of species to be simulated in the community.
}
  \item{family}{
A character string defining the error distribution and link function to be used when simulating the data. Only "probit" and "poisson" have been implemented so far.
}
  \item{paramX}{
A matrix of model parameters defining how each species (rows) is characterized by each descriptors (columns).
}
  \item{paramTr}{
A matrix of model parameters defining how each descriptors (rows) characterizes species traits (columns).
}
  \item{meansparamX}{
Vector of means used to generate the model parameters. (See details)
}
  \item{varX}{
Covariance matrix used to generate the model parameters. (See details)
}
  \item{latent}{
A list of matrices where each of set defines a random effect on the sites (rows) characterized by different latent variables (columns).
}
  \item{paramLatent}{
A list of matrices where each of set defines the parameters of a random effect on the sites (rows) characterized by different latent variables (columns).
}
  \item{spCor}{
A species by species correlation matrix defining the correlations among species in the community.
}
}
\details{
This function can be used to simulate a community from randomly proposed or fixed parameters. 

The descriptors in \code{X} are used without any modifications (or additions) to simulate the species community. As such, a column of 1 should be included in \code{X} for the model used to construct the community to include an intercept.

The values in \code{meansparamX} and \code{varX} are used as parameter of a multivariate normal distribution to generate the model's parameters (\code{\link[MASS]{mvrnorm}} is used in the function). When \code{paramX} is set to \code{NULL}, the \code{meansparamX} and the \code{varX} will be randomly generated if they are also set to \code{NULL}. When values are given to \code{paramX} the values of the \code{meansparamX} and the \code{varX} are not used and if either is set to \code{NULL}, no data will be generated for either set of parameter. When generated, the values for \code{meansparamX} are randomly sampled from \code{\link{rnorm}} and the values for \code{varX} are randomly sampled from an inverse (\code{solve})  Wishart distribution (\code{\link{rWishart}}).

Note that \code{meansparamX} can be calculated directly by \code{Tr} with \code{paramTr}. As such, \code{meansparamX} is made available as an extension. If it is given to the function but \code{paramTr} is not, than \code{paramTr} is calculated from \code{meansparamX} and vice versa. If both \code{meansparamX} and \code{paramTr} are given to the function and there is a mismatch in the parameters calculated, \code{meansparamX} will take precedent.

Currently, this function simulates two types of community:

	\tabular{ll}{
		\code{probit} \tab For presence/absence of species using a probit link.
		\code{poisson} \tab For count data of species using a log link.
	}
}
\value{
  The functions \code{communitySimulH} and \code{communitySimulHT} return an object of class \code{communitySimul} with the following components: 
  \item{data}{an object of class HMSCdata}
  \item{param}{an object of class HMSCparam}
  \item{probMat}{A matrix that defines the occurrence probability of each species for each site.}
}
\author{
F. Guillaume Blanchet
}

\examples{
### Construct some descriptors
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
### Simulate community from random parameters
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
desc<-cbind(1,scale(1:50),scale(1:50)^2)
traits <- rbind(1,1:20/5)
random1<-as.factor(1:50)
random2<-as.factor(rep(letters[1:2],each=25))
randEff<-data.frame(rand1=random1,rand2=random2)

### Simulate presence-absence community data
commDescTraitRandEffProbit<-communitySimul(X=desc,Tr=traits,
							Random=randEff,nsp=20)

### Simulate count community data (following poisson distribution)
commDescTraitRandEffPoisson<-communitySimul(X=desc,Tr=traits,
							 nsp=20,family="poisson")

#--------------------
### Fixing parameters
#--------------------
## ProbMat does not change
commParamX<-communitySimul(X=desc,
				paramX=commDescTraitRandEffProbit$param$paramX) 

# ProbMat changes
commParamTr<-communitySimul(X=desc,Tr=traits,
				paramTr=commDescTraitRandEffProbit$param$paramTr)

# ProbMat changes
commParamLatent<-communitySimul(X=desc,Random=randEff,
					paramLatent=commDescTraitRandEffProbit$param$paramLatent,
					latent=commDescTraitRandEffProbit$param$latent)

# ProbMat changes
speciesCor<-cov2cor(solve(rWishart(1,22,diag(20))[,,1]))
commParamSpCor<-communitySimul(X=desc,spCor=speciesCor)

}

\keyword{ datagen }
\keyword{ htest }
