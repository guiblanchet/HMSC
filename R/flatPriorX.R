#' @title Define flat priors 
#'
#' @description These functions define flat prior for the different models that can be constructed within the HMSC framework, they are meant to be used internally.
#'
#' @param data An object of the class \code{HMSCdata}.
#' @param meansParamX A vector of parameters that defines the prior information associated to parameters of the same name. 
#' @param varMeansParamX Symmetric covariance matrix. Prior for meansParamX defining how the sampling of the parameters vary.
#' @param varXDf Numeric. Prior for varX. It is the number of degrees of freedom of the inverse Wishart distribution from which varX is sampled.
#' @param varXScaleMat Symmetric covariance matrix. Prior for varX. It is the scale matrix of the inverse Wishart distribution from which varX is sampled.
#' @param paramTr Matrix of prior defining how descriptors (rows) characterizes traits (columns).
#' @param varTr Symmetric covariance matrix. Prior for paramTr defining how the sampling of the parameters vary.
#' @param paramPhylo Prior information associated to the parameter of the same name.
#' @param residVar Prior information associated to the parameter of the same name.
#' @param shrinkLocal Numeric. Hyperparameter of a gamma distribution defining the local shrinkage for latent variables. 
#' @param shrinkOverall A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the overall shrinkage of latent variable is handled. 
#' @param shrinkSpeed A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the speed of the shrinkage of latent variable is handled. 
#' @param shrinkLocalAuto List. Each level of the list includes a matrix of values describing how much shrinkage should be imposed on each parameter of the autocorrelated latent (unsampled) variables (\code{paramLatentAuto}). There are as many level in the list as there are columns in \code{Auto}. The size of each matrix is the same as the size of \code{paramLatentAuto}.
#' @param shrinkOverallAuto A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the overall shrinkage of the autocorrelated latent variable is handled. 
#' @param shrinkSpeedAuto A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the speed of the shrinkage of the autocorrelated latent variable is handled. 
#' @param paramAutoWeight Prior. Matrix with a single column defining the importance each value in \code{priorParamAutoDist} may have. The size of \code{priorParamAutoWeight} needs to be equal to the size of \code{priorParamAutoGrid}.
#' @param paramAutoDist Prior. Matrix with a single column defining the potential values given in \code{paramAuto} can take. The size of \code{priorParamAutoGrid} needs to be equal to the size of \code{priorParamAutoWeight}.
#' @param family Character string defining the type of generalized regression model for which priors need to be defined.
#' @param nparamX Numeric. Number of covariates.
#' @param nTr Numeric. Number of traits.
#'
#' @author F. Guillaume Blanchet
#'
#' @export
flatPriorX <- 
function(varXDf=NULL,varXScaleMat=NULL,meansParamX=NULL,varMeansParamX=NULL,nparamX){
	### Prior for varX
	if(is.null(varXDf)){
		varXDf<-nparamX+1
	}
	if(is.null(varXScaleMat)){
		varXScaleMat<-diag(nparamX)
	}
	
	### Prior for meansParamX
	if(is.null(meansParamX)){
		meansParamX<-as.matrix(rep(0,nparamX))
	}

	if(is.null(varMeansParamX)){
		varMeansParamX<-diag(nparamX)
	}
	
	### List of priors
	priors<-list(meansParamX=meansParamX,
				 varMeansParamX=varMeansParamX,
				 varXDf=varXDf,
				 varXScaleMat=varXScaleMat)
	
	return(priors)
}
