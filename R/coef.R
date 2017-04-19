#' @title Extract Model Coefficients for HMSC objects
#'
#' @description This function extract the average model coefficients from an \code{HMSC} object
#'
#' @param object Object of class hmsc.
#' @param \dots other arguments.
#'
#' @details
#'
#' This functions should be used when the MCMC run have converged properly for all parameters.
#'
#' The output of this function is an object of class \code{HMSCparam}.
#'
#' For models that include latent variables or autocorrelated latent variables, the latent variables and their associated parameters are also returned.
#'
#' @return
#'
#' An object of class \code{HMSCparam}.
#'
#' @author F. Guillaume Blanchet
#'
#' @seealso \code{\link{hmsc}}
#' @seealso \code{\link{predict.hmsc}}
#'
#' @examples
#' #================
#' ### Generate data
#' #================
#' desc <- cbind(scale(1:50), scale(1:50)^2)
#' nspecies <- 5
#' commDesc <- communitySimul(X = desc, nsp = nspecies)
#'
#' #=============
#' ### Formatting
#' #=============
#' ### Format data
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 2000, nburn = 1000, thin = 1)
#'
#' #==============================
#' ### Check for model convergence
#' #==============================
#' plot(as.mcmc(modelDesc,parameters = "paramX"))
#' plot(as.mcmc(modelDesc,parameters = "meansParamX"))
#' plot(as.mcmc(modelDesc,parameters = "varX"))
#'
#' #===============================
#' ### Extract average coefficients
#' #===============================
#' modelCoef <- coef(modelDesc)
#'
#' @keywords model
#' @keywords regression
#' @export
coef.hmsc <-
function(object,...){
#### F. Guillaume Blanchet
##########################################################################################
	### Build result object
	res <- vector("list",length=length(object$results$estimation))
	names(res) <- names(object$results$estimation)

	### Basic objects
	nAuto<-length(object$data$Auto)
	nRandom<-ncol(object$data$Random)
	nsp<-ncol(object$object$data$Y)

	### Number of iterations
	if(!is.null(nAuto)){
		niter<-length(object$results$estimation$paramLatentAuto)
	}
	if(!is.null(nRandom)){
		niter<-length(object$results$estimation$paramLatent)
	}
	if(!is.null(object$results$estimation$paramX)){
		niter<-dim(object$results$estimation$paramX)[3]
	}

	### paramX
	if(any(names(res)=="paramX")){
		res$paramX <- apply(object$results$estimation$paramX, 1:2, mean)
	}

	### meansParamX
	if(any(names(res)=="meansParamX")){
		res$meansParamX <- colMeans(object$results$estimation$meansParamX)
	}

	### paramTr
	if(any(names(res)=="paramTr")){
		res$paramTr <- apply(object$results$estimation$paramTr, 1:2, mean)
	}

	### varX
	if(any(names(res)=="varX")){
		res$varX <- apply(object$results$estimation$varX, 1:2, mean)
	}

	### paramPhylo
	if(any(names(res)=="paramPhylo")){
		res$paramPhylo <- mean(object$results$estimation$paramPhylo)
	}

	### paramLatent
	if(any(names(res)=="paramLatent")){
		dims <- sapply(object$results$estimation$paramLatent,dim)

		paramLatent <- array(0,dim=c(max(dims[1,]),max(dims[2,]),niter))

		for(i in 1:niter){
			paramLatent[1:dims[1,i],1:dims[2,i],] <- object$results$estimation$paramLatent[[i,1]]
		}
		res$paramLatent <- apply(paramLatent, 1:2, mean)
	}

	### latent
	if(any(names(res)=="latent")){
		dims <- sapply(object$results$estimation$latent,dim)

		latent <- array(0,dim=c(max(dims[1,]),max(dims[2,]),niter))

		for(i in 1:niter){
			latent[1:dims[1,i],1:dims[2,i],] <- object$results$estimation$latent[[i,1]]
		}
		res$latent <- apply(latent, 1:2, mean)
	}

	### paramLatentAuto
	if(any(names(res)=="paramLatentAuto")){
		dims <- sapply(object$results$estimation$paramLatentAuto,dim)

		paramLatentAuto <- array(0,dim=c(max(dims[1,]),max(dims[2,]),niter))

		for(i in 1:niter){
			paramLatentAuto[1:dims[1,i],1:dims[2,i],] <- object$results$estimation$paramLatentAuto[[i,1]]
		}
		res$paramLatentAuto <- apply(paramLatentAuto, 1:2, mean)
	}

	### latentAuto
	if(any(names(res)=="latentAuto")){
		dims <- sapply(object$results$estimation$latentAuto,dim)

		latentAuto <- array(0,dim=c(max(dims[1,]),max(dims[2,]),niter))

		for(i in 1:niter){
			latentAuto[1:dims[1,i],1:dims[2,i],] <- object$results$estimation$latentAuto[[i,1]]
		}
		res$latentAuto <- apply(latentAuto, 1:2, mean)
	}

	### paramAuto
	if(any(names(res)=="paramAuto")){
		paramAuto <- matrix(NA,nrow=niter,ncol=nrow(object$results$estimation$paramAuto[[1,1]]))

		for(i in 1:niter){
			paramAuto[i,] <- object$results$estimation$paramAuto[[i,1]]
		}

		res$paramAuto <- colMeans(paramAuto)
	}
	return(res)
}
