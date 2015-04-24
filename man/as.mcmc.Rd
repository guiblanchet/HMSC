\name{as.mcmc}
\alias{as.mcmc}
\title{
Convert MCMC parameters estimation to \code{\link{mcmc}[coda]} object 
}
\description{
Use the Markov Chain Monte Carlo results obtained from \code{\link{hmscProbit}} and convert them to an \code{\link{mcmc}[coda]} object.  
}
\usage{
as.mcmc(hmsc,parameters="paramX",burning=FALSE)
}
\arguments{
  \item{hmsc}{
Object of class hmsc.
}
  \item{parameters}{
A character string defining the parameters for which the convergions needs to be carried out.
}
  \item{burning}{
Logical. Whether the burning iterations should be include (\code{TRUE}) or not (\code{FALSE}).
}
}
\details{
The parameters estimated for the latent variables (\code{paramLatent}) cannot be directly converted to an object of class \code{mcmc} because the procedure implemented in \code{\link{hmscProbit}} allows for the number of latent variables, and thus of parameters associated to these latent variables, to change. For this reason, if \code{paramLatent} is used, it is the matrix defining the correlations among species calculated from  \code{paramLatent} for which an \code{\link{mcmc}[coda]} object is constructed. 

The \code{\link{mcmc}[coda]} object associated to \code{paramLatent} includes information only for the correlations among single pairs of species, that is, only the lower triangle of the correlation matrix is considered.

When there are multiple random effects estimated for the model, associated to each random effect there is a set of estimated latent parameters (\code{paramLatent}). For this reason, when the Markov Chain Monte Carlo results for \code{paramLatent} is converted with \code{as.mcmc}, parameters associated to each set of latent variables in \code{paramLatent} becomes an code{\link{mcmc}[coda]} object in a list. It is important to note that this "list" is \emph{NOT} an \code{mcmc.list}. 

}

\value{
An object of class \code{\link{mcmc}[coda]}.
}

\author{
F. Guillaume Blanchet
}
\seealso{
\code{\link{mcmc}[coda]}, \code{\link{hmscProbit}}
}
\examples{
#==================
### Simulating data
#==================
desc<-cbind(1,scale(1:50),scale(1:50)^2)
nspecies<-20
commDesc<-communitySimulProbit(X=desc,nsp=nspecies)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDesc$data$Y,X=desc,interceptX=FALSE,
					  interceptTr=FALSE)

#==============
### Build model
#==============
modelDesc<-hmscProbit(formdata,niter=200,nburn=100,thin=1,verbose=FALSE)

paramXmcmc<-as.mcmc(modelDesc,parameters="paramX")
}
\keyword{ IO }
