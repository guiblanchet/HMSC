\name{as.HMSC}
\alias{as.HMSCdata}
\alias{as.HMSCparam}
\alias{as.HMSCprior}
\title{
Format data, parameters, and priors to "HMSC" classes 
}
\description{
These functions are used to format the data, parameters, and priors to object of class "HMSC". 
}
\usage{
as.HMSCdata(Y,X,Tr=NULL,Random=NULL,scaleX=TRUE,scaleTr=TRUE,
			interceptX=TRUE,interceptTr=TRUE)
as.HMSCprior(HMSCdata,varXDf=NULL,varXScaleMat=NULL,paramTr=NULL,
			 varTr=NULL,shrinkOverall=NULL,shrinkSpeed=NULL,
			 shrinkLocal=NULL)
as.HMSCparam(HMSCdata,HMSCprior,paramX=NULL,varX=NULL,paramTr=NULL,
			 varTr=NULL,family=binomial(link = "probit"))
}
\arguments{
  \item{HMSCdata}{
An object of the class \code{HMSCdata}
}
  \item{HMSCprior}{
An object of the class \code{HMSCprior}. 
}
  \item{Y}{
A matrix or a data.frame where each column is a species.
}
  \item{X}{
A matrix or a data.frame where each column is a descriptor of the sites.
}
  \item{Tr}{
A matrix or a data.frame where each row is a trait characterizing all species.
}
  \item{Random}{
A factor defining a random effect on the samples or a data.frame characterizing multiple random effects.
}
  \item{scaleX}{
Logical. Whether the columns of X should be centred and divided by the standard deviation. Default is TRUE.
}
  \item{scaleTr}{
Logical. Whether the rows of Tr should be centred and divided by the standard deviation. Default is TRUE.
}
  \item{interceptX}{
Logical. Whether a column of 1s should be added to X. Default is TRUE.
}
  \item{interceptTr}{
Logical. Whether a row of 1s should be added to Tr. Default is TRUE.
}
  \item{paramX}{
A matrix or data.frame of model parameters defining how species (rows) are characterized by the descriptors (columns).
}
  \item{varX}{
Symmetric covariance matrix. Each dimension of this matrix should be equal to the number of explanatory variables.
}
  \item{paramTr}{
In function as.HMSCdata, it is a matrix of model parameters defining how descriptors (rows) characterizes traits (columns).

In function as.HMSCprior, it is a matrix of prior defining how descriptors (rows) characterizes traits (columns).
}
  \item{varTr}{
In function as.HMSCdata, it is a symmetric covariance matrix. Each dimension of this matrix should be equal to the number of traits.

In function as.HMSCprior, it is a symmetric covariance matrix of prior. Each dimension of this matrix should be equal to the number of traits.
}
  \item{family}{
Scalar defining the within group variance in \code{Random}. (See details)
}
  \item{varXDf}{
Numeric. Prior for varX. It is the number of degrees of freedom of the inverse Wishart distribution from which varX is sampled.
}
  \item{varXScaleMat}{
Symmetric covariance matrix. Prior for varX. It is the scale matrix of the inverse Wishart distribution from which varX is sampled.
}
  \item{shrinkOverall}{
A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the overall shrinkage of latent variable is handled (see details). 
}
  \item{shrinkSpeed}{
A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the speed of the shrinkage of latent variable is handled (see details). 
}
  \item{shrinkLocal}{
Numeric. Hyperparameter of a gamma distribution defining the local shrinkage for latent variables. 
}
}
\details{
Functions \code{as.HMSCparam} and \code{as.HMSCprior} do not need to be used when estimating a model using \code{\link{hmscProbit}}. In \code{\link{hmscProbit}}, if \code{as.HMSCparam} and \code{as.HMSCprior} are not defined, the function will define a set of default parameters and priors. \code{as.HMSCparam} and \code{as.HMSCprior} were designed so that some, not necessarily all, priors and parameters need to be defined. For example, if one defines only the prior (parameters) for \code{varTr}, all other priors (parameters) will be defined with default priors (parameter). 
 
In \code{as.HMSCparam}, the argument \code{family} is used to define the  \code{paramX} by default using a univariate generalized linear model for each species. Currently, only \code{binomial(link = "probit")} should be used because \code{\link{hmscProbit}} only estimate probit models. 

The choosing the values for \code{shrinkOverall} and \code{shrinkSpeed} is not straightforward. Bhattacharya and Dunson (2011) proposed to use a 2.1 and 1 for \code{shrinkOverall} and 3.1 and 1 for \code{shrinkSpeed}. However, for presence-absence data, it seems that these shrinkage values (especially the scale [first] parameter) are not large enough, which leads to a too many latent variables defined. It is the reason why we suggest to use larger scale values. As such, if they are not predefine, the scale value for \code{shrinkOverall} is 10 whereas the scale value for \code{shrinkSpeed} is 15. These two scale parameters are likely the only one that need to be chosen more carefully to ensure that the model estimation is accurate.
}

\value{
An object of the class \code{HMSCdata} is returned by \code{as.HMSCdata}.
An object of the class \code{HMSCparam} is returned by \code{as.HMSCparam}. 
An object of the class \code{HMSCprior} is returned by \code{as.HMSCprior}. 

All of these objects present the same output data as the input except for \code{as.HMSCparam} where \code{meansparamX} is also defined.
}
\references{
Bhattacharya, A. and Dunson, D.B. (2011) Sparse Bayesian infinite factor models. \emph{Biometrika} \strong{98}, 291--306.
}

\author{
F. Guillaume Blanchet
}
\note{
These functions will likely get more refined with time and with the new models implemented in the package. 
}

\examples{
#================
### Generate data
#================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
randomEff<-as.factor(rep(letters[1:5],each=10))
dataBase<-communitySimulProbit(desc,nsp=30,Random=randomEff)

#=============
### Formatting
#=============
### Format data
formdata<-as.HMSCdata(Y=dataBase$data$Y,X=desc,Random=randomEff,
					  interceptX=FALSE)

### Format priors
formpriors<-as.HMSCprior(formdata,shrinkOverall=c(50,1),
						 shrinkSpeed=c(20,1))

### Format parameters
formparam<-as.HMSCparam(formdata,formpriors,
						paramX=dataBase$param$paramX)
}
\keyword{manip}
\keyword{classes}
