#' @title Initiate a latent Y
#'
#' @description Initiate a latent set of response variables for the different model types
#'
#' @param data An object of the class \code{HMSCdata}.
#' @param param An object of the class \code{HMSCparam}.
#' @param family A character string defining how the latent response variables needs to be initiated.
#'
#' @details
#'
#' \code{HMSCdata} and \code{HMSCparam}, the proper set of starting parameters are generated.
#' Since the function is meant to be used by expert users, there are no predefined setting for any of the argument of the function.
#'
#' @return
#'
#'  A matrix or a list of objects that are needed to construct the model of interest.
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @keywords datagen
#' @export
iniYlatent <-
function(data,param,family){

	### A few basic objects
	if(!is.null(data$X)){
		nparamX<-ncol(data$X)
	}
	nsp<-ncol(data$Y)
	nsite<-nrow(data$Y)

	#==========================================
	### Initiate a latent Y for gaussian models
	#==========================================
	if(family=="gaussian"){
		### Transform Y into a matrix
		Y<-as.matrix(data$Y)

		### Initial latent Y
		Ylatent<-Y
	}

	#=======================================
	### Initiate a latent Y for other models
	#=======================================
	if(family!="gaussian"){
		### Transform X and Y into a matrix
		Y<-as.matrix(data$Y)
		if(!is.null(data$X)){
			X<-as.matrix(data$X)
		}

		#---------------------
		### Calculate EstModel
		#---------------------
		if(is.null(data$X)){
			if(is.null(data$Random)){
				if(is.null(data$Auto)){
					EstModel<-matrix(0,nsite,nsp)
				}else{
					### Only Auto
					EstModel<-matrix(0,nsite,nsp)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
				}
			}else{
				if(is.null(data$Auto)){
					### Only Random
					EstModel<-matrix(0,nsite,nsp)
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}else{
					### Random and Auto
					EstModel<-matrix(0,nsite,nsp)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}
			}
		}else{
			if(is.null(data$Random)){
				if(is.null(data$Auto)){
					### Only X
					EstModel<-tcrossprod(X,param$param$paramX)
				}else{
					### X and Auto
					EstModel<-tcrossprod(X,param$param$paramX)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
				}
			}else{
				if(is.null(data$Auto)){
					### X and Random
					EstModel<-tcrossprod(X,param$param$paramX)
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}else{
					### X, Random and Auto
					EstModel<-tcrossprod(X,param$param$paramX)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}
			}
		}

		if(family == "probit"){
			Y0<-Y==0
			Y1<-Y==1

			Ylatent<-sampleYlatentProbit(Y0,Y1,matrix(0,nrow=nsite,ncol=nsp),EstModel,param$param$residVar,nsp,nsite)
		}

		if(family == "poisson" | family == "overPoisson"){
			Ylatent<-sampleYlatentPoisson(Y,matrix(0,nrow=nsite,ncol=nsp),EstModel,param$param$residVar,nsp,nsite)
		}

	}

	class(Ylatent)<-"HMSCYlatent"
	return(Ylatent)
}
