\name{hmsc}
\alias{hmsc}
\title{
Hierarchical modelling of species community
}
\description{
Hierarchical modelling of species community using species occurrence (presence-absence) data. 
}
\usage{
hmsc(data,param=NULL,priors=NULL,family="probit",niter=1000,nburn=500,
			thin=1,verbose=TRUE)
}
\arguments{
  \item{data}{
Object of class HMSCdata.
}
  \item{param}{
Object of class HMSCparam.
}
  \item{priors}{
Object of class HMSCprior.
}
  \item{family}{
A character string defining the error distribution and link function to be used when performing the analysis. Only "probit" and "poisson" have been implemented so far.
}
  \item{niter}{
Numeric. A positive integer defining the number of iterations to be carried out. Default is 1000.
}
  \item{nburn}{
Numeric. A positive integer defining the number of iterations to be used in the burning (first) phase of the algorithm. Default is 500.
}
  \item{thin}{
Numeric. A positive integer defining thinning. Default is 1 (see details).
}
  \item{verbose}{
Logical. Whether a counter printing on the screen how many iterations have been performed so far is used or not. Default is TRUE.
}
}
\details{
If \code{param} and \code{priors} are not defined in the function, default parameters and priors will be defined.  

When thinning, the value given to \code{thin} means that all but the \emph{k}^{th} values will be considered.  

This function makes it possible to build models that consider, habitat characteristics, traits and random effects (one or multiple) on the intercept. The way to define how these models relies on the data included in \code{data}.
}

\value{
The output of this function varies according the model built. However, they follow the same general structure. The output is always divided in a \code{burning} and an \code{estimation} part which contain the following parts: 
	\item{\code{paramX}}{A 3-dimensional array where, for each species (first dimension), the regression parameters of X (second dimension) are stored for all iterations (third dimension).}
	\item{\code{paramTr}}{A 3-dimensional array where, for each regression parameter (first dimension), the species traits (second dimension) are stored for all iterations (third dimension).}
	\item{\code{varX}}{A 3-dimensional array where covariance matrices measuring how the species vary in their responses to the descriptors (first and second dimension) are stored for all iterations (third dimension).}
	\item{\code{latent}}{A list of lists (one for each random effect) that includes matrices of different dimensions. Each set of lists corresponds to an iteration, the row of the matrices corresponds to samples and the columns to latent variables }
	\item{\code{paramLatent}}{A list of lists (one for each random effect) that includes matrices of different dimensions. Each set of lists corresponds to an iteration, the row of the matrices corresponds to species and the columns to latent variables.}
}
\references{
Bhattacharya, A. and Dunson, D.B. (2011) Sparse Bayesian infinite factor models. \emph{Biometrika} \strong{98}, 291--306.

Ovaskainen, O. and Soininen, J. (2011) Making more out of sparse data: hierarchical modeling of species communities. \emph{Ecology} \strong{92}, 289-295.
}
\author{
F. Guillaume Blanchet
}
\seealso{
\code{\link{as.HMSCdata}}, \code{\link{as.HMSCparam}},\code{\link{as.HMSCprior}}, \code{\link{as.mcmc}}

}
\examples{
#################################################################
### Examples for presence-absence data using Probit link function
#################################################################
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Only site characteristics
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
nspecies<-20
commDesc<-communitySimul(X=desc,nsp=nspecies)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDesc$data$Y,X=desc,interceptX=FALSE,
					  interceptTr=FALSE)

#==============
### Build model
#==============
modelDesc<-hmsc(formdata,niter=200,nburn=100,thin=1,verbose=FALSE)

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Only traits characteristics
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
traits <- rbind(1,1:20/5)
nspecies<-20

commTrait<-communitySimul(X=desc,Tr=traits,nsp=nspecies)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commTrait$data$Y,X=desc,Tr=traits,
					  interceptX=FALSE,interceptTr=FALSE)

#==============
### Build model
#==============
modelTrait<-hmsc(formdata,niter=200,nburn=100,thin=1,verbose=FALSE)

#_+_+_+_+_+_+_+_+_+_+_+_+_
### Only one random effect
#_+_+_+_+_+_+_+_+_+_+_+_+_
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
randEff<-as.factor(rep(letters[1:5],each=10))
nspecies<-20

commRandEff1<-communitySimul(X=desc,Random=randEff,nsp=nspecies)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commRandEff1$data$Y,X=desc,Random=randEff,
					  interceptX=FALSE)

#==============
### Build model
#==============
modelRandEff1<-hmsc(formdata,niter=200,nburn=100,thin=1,
						  verbose=FALSE)

#_+_+_+_+_+_+_+_+_+_+_+_+_
### Only two random effect
#_+_+_+_+_+_+_+_+_+_+_+_+_
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
random1<-as.factor(1:50)
random2<-as.factor(rep(letters[1:2],each=25))
randEff<-data.frame(rand1=random1,rand2=random2)
nspecies<-20

commRandEff2<-communitySimul(X=desc,Random=randEff,nsp=nspecies)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commRandEff2$data$Y,X=desc,Random=randEff,
					  interceptX=FALSE)

#==============
### Build model
#==============
modelRandEff2<-hmsc(formdata,niter=200,nburn=100,thin=1,
						  verbose=FALSE)

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Site characteristics and traits
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
traits <- rbind(1,1:20/5)
commDescTrait<-communitySimul(X=desc,Tr=traits)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDescTrait$data$Y,X=desc,Tr=traits,
					  interceptX=FALSE,interceptTr=FALSE)

#==============
### Build model
#==============
modelDescTrait<-hmsc(formdata,niter=200,nburn=100,thin=1,
						   verbose=FALSE)


#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Site characteristics and 2 random effects
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
random1<-as.factor(1:50)
random2<-as.factor(rep(letters[1:2],each=25))
randEff<-data.frame(rand1=random1,rand2=random2)
commDescRandEff<-communitySimul(X=desc,Random=randEff,nsp=20)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDescRandEff$data$Y,X=desc,Random=randEff,
					  interceptX=FALSE)

#==============
### Build model
#==============
modelDescRandEff<-hmsc(formdata,niter=200,nburn=100,thin=1,
							 verbose=FALSE)

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Site characteristics, traits and 2 random effects
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
traits <- rbind(1,1:20/5)
random1<-as.factor(1:50)
random2<-as.factor(rep(letters[1:2],each=25))
randEff<-data.frame(rand1=random1,rand2=random2)
commDescTraitRandEff<-communitySimul(X=desc,Tr=traits,
								Random=randEff,nsp=20)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDescTraitRandEff$data$Y,X=desc,Tr=traits,
					  Random=randEff,interceptX=FALSE)

#==============
### Build model
#==============
modelDescTraitRandEff<-hmsc(formdata,niter=200,nburn=100,thin=1,
								  verbose=FALSE)

##################################################################
### Examples for presence-absence data using Poisson link function
##################################################################
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Only site characteristics
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
nspecies<-20
commDesc<-communitySimul(X=desc,nsp=nspecies,family="poisson")

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDesc$data$Y,X=desc,interceptX=FALSE,
					  interceptTr=FALSE)

#==============
### Build model
#==============
modelDesc<-hmsc(formdata,niter=200,nburn=100,thin=1,verbose=FALSE,family="poisson")

#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
### Only traits characteristics
#_+_+_+_+_+_+_+_+_+_+_+_+_+_+_+
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
traits <- rbind(1,1:20/5)
nspecies<-20

commTrait<-communitySimul(X=desc,Tr=traits,nsp=nspecies,family="poisson")

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commTrait$data$Y,X=desc,Tr=traits,
					  interceptX=FALSE,interceptTr=FALSE)

#==============
### Build model
#==============
modelTrait<-hmsc(formdata,niter=200,nburn=100,thin=1,verbose=FALSE,family="poisson")

#_+_+_+_+_+_+_+_+_+_+_+_+_
### Only one random effect
#_+_+_+_+_+_+_+_+_+_+_+_+_
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
randEff<-as.factor(rep(letters[1:5],each=10))
nspecies<-20

commRandEff1<-communitySimul(X=desc,Random=randEff,nsp=nspecies,family="poisson")

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commRandEff1$data$Y,X=desc,Random=randEff,
					  interceptX=FALSE)

#==============
### Build model
#==============
modelRandEff1<-hmsc(formdata,niter=200,nburn=100,thin=1,
						  verbose=FALSE,family="poisson")

#_+_+_+_+_+_+_+_+_+_+_+_+_
### Only two random effect
#_+_+_+_+_+_+_+_+_+_+_+_+_
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
random1<-as.factor(1:50)
random2<-as.factor(rep(letters[1:2],each=25))
randEff<-data.frame(rand1=random1,rand2=random2)
nspecies<-20

commRandEff2<-communitySimul(X=desc,Random=randEff,nsp=nspecies,family="poisson")

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commRandEff2$data$Y,X=desc,Random=randEff,
					  interceptX=FALSE)

#==============
### Build model
#==============
modelRandEff2<-hmsc(formdata,niter=200,nburn=100,thin=1,
						  verbose=FALSE,family="poisson")
}
\keyword{ multivariate }
\keyword{ models }
