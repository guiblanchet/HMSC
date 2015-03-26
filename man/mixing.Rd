\name{mixing}
\alias{mixing}
\title{
Mixing plot
}
\description{
This function draws mixing plots for the parameters estimated.
}
\usage{
mixing(x, toplot = "paramX", param = c(1, 1), add = FALSE, ...)
}
\arguments{
  \item{x}{
Object of class hmsc
}
  \item{toplot}{
A character string defining whether the mixing plot should be drawn for parameter in "paramX", "paramTr", or "varX". Any unambiguous variation of wording used in this argument is accepted.
}
  \item{param}{
A vector of one two integer defining which parameter should the mixing plot be drawn. 
}
  \item{add}{
Logical. Whether a new plot (\code{FALSE}) should be drawn or the mixing plot should be drawn within another plot.
}
  \item{\dots}{
Graphical arguments passed to this function.
}
}
\author{
F. Guillaume Blanchet
}
\examples{
desc<-cbind(1,scale(1:50),scale(1:50)^2)
traits <- rbind(1,1:20/5)
commDescTrait<-communitySimulProbit(X=desc,Tr=traits)

#=================================
### Formatting data and parameters
#=================================
formdata<-as.HMSCdata(Y=commDescTrait$data$Y,X=desc,Tr=traits,
					  interceptX=FALSE,interceptTr=FALSE)

#==============
### Build model
#==============
modelDescTrait<-hmscProbit(formdata,niter=200,nburn=100,thin=1,
							verbose=FALSE)

#=========
### Mixing
#=========
#---------
### paramX
#---------
par(mfrow=c(5,4),mar=c(1,3,2,1))
for(i in 1:20){
	mixing(modelDescTrait,toplot="paramX",param=c(i,1),
		   ylim=range(modelDescTrait$estimation$paramX[i,,]),
		   xlab="Sites",ylab="paramX",las=1,xaxt="n")
	title(paste("Species",i))
	mixing(modelDescTrait,toplot="paramX",param=c(i,2),add=TRUE,
		   col="blue")
	mixing(modelDescTrait,toplot="paramX",param=c(i,3),add=TRUE,
		   col="red")
}

#-------
### varX
#-------
par(mfrow=c(3,3),mar=c(1,3,3,1))
for(i in 1:3){
	for(j in 1:3){
		mixing(modelDescTrait,toplot="varX",param=c(i,j),col="red",
			   xlab="Sites",ylab="Sigma",las=1)
		title(paste("Community varX",i,"-",j))
	}
}

#----------
### ParamTr
#----------
par(mfrow=c(3,2),mar=c(1,3,3,1))
for(i in 1:3){
	for(j in 1:2){
		mixing(modelDescTrait,toplot="paramTr",param=c(i,j),col="red",
			   xlab="Sites",ylab="Sigma",las=1)
		title(paste("Community varX",i,"-",j))
	}
}
}
\keyword{ hplot }
\keyword{ aplot }
