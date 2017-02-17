#' @title Initiate parameters
#'
#' @description These functions are used to initiate the different sets of parameters used to construct different type of models within the HMSC framework, they are meant to be used internally.
#'
#' @param data An object of the class \code{HMSCdata}.
#' @param priors An object of the class \code{HMSCprior}. 
#' @param paramX A matrix or data.frame of model parameters defining how species (rows) are characterized by the descriptors (columns).
#' @param meansParamX A vector of parameters defining how an average species reacts to each descriptor in X. 
#' @param varX Symmetric covariance matrix. Each dimension of this matrix should be equal to the number of explanatory variables.
#' @param paramTr Matrix of model parameters defining how descriptors (rows) characterizes traits (columns). 
#' @param paramPhylo Numeric. Defines the importance of phylogeny in structuring the community. 
#' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model. 
#' @param latent List. Each level of the list includes a matrix of latent (unsampled) variables. There are as many level in the list as there are columns in \code{Random}. The number of rows in each matrix corresponds to the number of levels in each factor in \code{Random}. 
#' @param paramLatent List. Each level of the list includes a matrix of parameters for the latent (unsampled) variables. There are as many level in the list as there are columns in \code{Random}. Each matrix has as many rows as there are species in \code{Y} and as many columns as there are latent variables in \code{latent}.
#' @param shrinkLocal List. Each level of the list includes a matrix of values describing how much shrinkage should be imposed on each parameter of the latent (unsampled) variables (\code{paramLatent}). There are as many level in the list as there are columns in \code{Random}. The size of each matrix is the same as the size of \code{paramLatent}.
#' @param paramShrinkGlobal List. Each level of the list includes a vector of values describing how much shrinkage should be imposed on each latent variables (\code{latent}) as a whole. There are as many level in the list as there are columns in \code{Random}. The number of elements in each vector equals the number of variables in \code{latent}.
#' @param paramAuto Vector of numerical values defining the autocorrelated patterns for each factor in \code{Auto}. These values can range from 0 to the largest distance between samples in the autocorrelated level. 
#' @param latentAuto List. Each level of the list includes a matrix of autocorrelated latent (unsampled) variables. There are as many level in the list as there are columns in \code{Auto}. The number of rows in each matrix corresponds to the number of levels in each autocorrelated factor in \code{Auto}. 
#' @param paramLatentAuto List. Each level of the list includes a matrix of parameters for the autocorrelated latent (unsampled) variables. There are as many level in the list as there are columns in \code{Auto}. Each matrix has as many rows as there are species in \code{Y} and as many columns as there are autocorrelated latent variables in \code{latentAuto}.
#' @param shrinkLocalAuto List. Each level of the list includes a matrix of values describing how much shrinkage should be imposed on each parameter of the autocorrelated latent (unsampled) variables (\code{paramLatentAuto}). There are as many level in the list as there are columns in \code{Auto}. The size of each matrix is the same as the size of \code{paramLatentAuto}.
#' @param paramShrinkGlobalAuto List. Each level of the list includes a vector of values describing how much shrinkage should be imposed on each autocorrelated latent variables (\code{latent}) as a whole. There are as many level in the list as there are columns in \code{Auto}. The number of elements in each vector equals the number of variables in \code{latentAuto}.
#'
#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats rWishart
#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats binomial
#' @importFrom stats lm
#'
#' @author F. Guillaume Blanchet
#'
#' @export
iniParamX <-
function(data,priors,paramX=NULL,meansParamX=NULL,varX=NULL){
	
	### Number of X variables
	nparamX<-ncol(data$X)
	
	### Number of species
	nsp<-ncol(data$Y)
	
	#------------------
	### Initiate paramX
	#------------------
	if(nrow(data$Y)!=0){
		
		if(is.null(paramX)){
			paramX<-matrix(NA,nrow=nsp,ncol=nparamX)
			
			options(warn=-1)
			if(attributes(priors)$distr=="probit"){
				for(i in 1:nsp){
					paramX[i,]<-coef(glm(data$Y[,i]~-1+.,data=as.data.frame(data$X),family=binomial(link = "probit")))
				}
			}
			if(attributes(priors)$distr=="gaussian"){
				for(i in 1:nsp){
					paramX[i,]<-coef(lm(data$Y[,i]~-1+.,data=as.data.frame(data$X)))
				}
			}
			options(warn=0)
			
			### Correct for extreme values
			paramXtoCorrPos<-which(paramX>4,arr.ind=TRUE)
			paramXtoCorrNeg<-which(paramX< -4,arr.ind=TRUE)
			
			if(length(paramXtoCorrPos)>0){
				paramX[paramXtoCorrPos]<- 4
			}
			if(length(paramXtoCorrNeg)>0){
				paramX[paramXtoCorrNeg]<- -4
			}
		}
	}else{
		paramX<-matrix(rnorm(nsp*nparamX),nrow=nsp,ncol=nparamX)
	}
	#-----------------------
	### Initiate meansParamX
	#-----------------------
	if(is.null(meansParamX)){
		meansParamX<-1/(nsp+1)*(colSums(paramX)+priors$param$meansParamX)
	}
	#----------------
	### Initiate varX
	#----------------
	if(is.null(varX)){
		varX<-cov(paramX)+diag(0.1,nparamX)
		
		if(any(is.na(varX))){
			precX <- rWishart(1,paramX+1,diag(nparamX))[,,1]
			varX <- solve(precX)
		}else{
			precX<-solve(varX)
		}
	}else{
		precX<-solve(varX)
	}
	
	
	### List of parameters
	param<-list(paramX=paramX,
				varX=varX,
				precX=precX,
				meansParamX=meansParamX)
	
	### Name objects
	rownames(param$paramX)<-colnames(data$Y)
	colnames(param$paramX)<-colnames(data$X)
	
	rownames(param$varX)<-colnames(data$X)
	colnames(param$varX)<-colnames(data$X)
	
	rownames(param$precX)<-colnames(data$X)
	colnames(param$precX)<-colnames(data$X)
	
	rownames(param$meansParamX)<-colnames(data$X)
	
	return(param)
}
