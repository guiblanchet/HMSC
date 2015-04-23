\name{corRandomEff}
\alias{corRandomEff}
\title{
Construct a correlation matrix from the parameters of the latent variables (\code{paramLatent})
}
\description{
Construct a correlation matrix from the parameters of the latent variables (\code{paramLatent})
}
\usage{
corRandomEff(hmsc,burning=FALSE)
}
\arguments{
  \item{hmsc}{
Object of class hmsc.
}
  \item{burning}{
Logical. Whether the burning iterations should be include (\code{TRUE}) or not (\code{FALSE}).
}
}

\details{
The output of this function is an array  object 
}

\value{
An object of class \code{corRandomEff}. This object is an array of four dimensions. The correlation between two species (first and second dimension) can be found for a particular iteration (third dimension) associated to a specific random effect.
}

\author{
F. Guillaume Blanchet
}
\seealso{
\code{\link{hmscProbit}}
}
\examples{
#==================
### Simulating data
#==================
desc<-matrix(1,nrow=50,ncol=1)
random1<-as.factor(1:50)
random2<-as.factor(rep(letters[1:2],each=25))
randEff<-data.frame(rand1=random1,rand2=random2)
nspecies<-10

commRandEff2<-communitySimulProbit(X=desc,Random=randEff,nsp=nspecies)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commRandEff2$data$Y,X=desc,Random=randEff,
					  interceptX=FALSE)

#==============
### Build model
#==============
modelRandEff2<-hmscProbit(formdata,niter=200,nburn=100,thin=1,
						  verbose=FALSE)

#===============================				  
### Construct correlation matrix
#===============================				  
corRes<-corRandomEff(modelRandEff2)

#=================================				  
### Draw chord diagram with corRes
#=================================				  
library(circlize)

chordDiagram(corRes[,,1,1],symmetric=TRUE)

}
\keyword{ IO }
