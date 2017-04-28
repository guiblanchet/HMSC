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
		### Number of random effect
		nRandom <- ncol(object$results$estimation$paramLatent)

		### Number of latent variables for each random effect
		dims <- array(NA, dim=c(2,niter, nRandom))

		for(i in 1:nRandom){
			dims[,,i] <- sapply(object$results$estimation$paramLatent[,i],dim)
		}

		### Store paramLatent for each random effect
		paramLatent <- vector("list",length=nRandom)
		for(i in 1:nRandom){
			paramLatent[[i]] <- array(0,dim=c(max(dims[1,,i]),max(dims[2,,i]),niter))
		}

		for(i in 1:nRandom){
			for(j in 1:niter){
				paramLatent[[i]][1:dims[1,j,i],1:dims[2,j,i],j] <- object$results$estimation$paramLatent[[j,i]]
			}
		}

		### Calculate the average coefficients
		paramLatentMean <- vector("list",length=nRandom)
		for(i in 1:nRandom){
			paramLatentMean[[i]]<-apply(paramLatent[[i]], 1:2, mean)
		}

		### names each part of the object
		for(i in 1:nRandom){
			rownames(paramLatentMean[[i]]) <- rownames(object$results$estimation$paramLatent[[1,i]])
			colnames(paramLatentMean[[i]]) <- colnames(object$results$estimation$paramLatent[[which.max(dims[2,,i]),i]])
		}

		paramLatentMean<-matrix(paramLatentMean,ncol=nRandom)
		colnames(paramLatentMean)<-colnames(object$results$estimation$paramLatent)

		res$paramLatent <- paramLatentMean
	}

	### latent
	if(any(names(res)=="latent")){
		### Number of random effect
		nRandom <- ncol(object$results$estimation$latent)

		### Number of latent variables for each random effect
		dims <- array(NA, dim=c(2,niter, nRandom))

		for(i in 1:nRandom){
			dims[,,i] <- sapply(object$results$estimation$latent[,i],dim)
		}

		### Store paramLatent for each random effect
		latent <- vector("list",length=nRandom)
		for(i in 1:nRandom){
			latent[[i]] <- array(0,dim=c(max(dims[1,,i]),max(dims[2,,i]),niter))
		}

		for(i in 1:nRandom){
			for(j in 1:niter){
				latent[[i]][1:dims[1,j,i],1:dims[2,j,i],j] <- object$results$estimation$latent[[j,i]]
			}
		}

		### Calculate the average coefficients
		latentMean <- vector("list",length=nRandom)
		for(i in 1:nRandom){
			latentMean[[i]]<-apply(latent[[i]], 1:2, mean)
		}

		### names each part of the object
		for(i in 1:nRandom){
			rownames(latentMean[[i]]) <- rownames(object$results$estimation$latent[[1,i]])
			colnames(latentMean[[i]]) <- colnames(object$results$estimation$latent[[which.max(dims[2,,i]),i]])
		}

		latentMean<-matrix(latentMean,ncol=nRandom)
		colnames(latentMean)<-colnames(object$results$estimation$latent)

		res$latent <- latentMean
	}

	### paramLatentAuto
	if(any(names(res)=="paramLatentAuto")){
		### Number of autocorrelated random effect
		nAuto <- ncol(object$results$estimation$paramLatentAuto)

		### Number of autocorrelated latent variables for each autocorrelated random effect
		dims <- array(NA, dim=c(2,niter, nAuto))

		for(i in 1:nAuto){
			dims[,,i] <- sapply(object$results$estimation$paramLatentAuto[,i],dim)
		}

		### Store paramLatentAuto for each autocorrelated random effect
		paramLatentAuto <- vector("list",length=nAuto)
		for(i in 1:nAuto){
			paramLatentAuto[[i]] <- array(0,dim=c(max(dims[1,,i]),max(dims[2,,i]),niter))
		}

		for(i in 1:nAuto){
			for(j in 1:niter){
				paramLatentAuto[[i]][1:dims[1,j,i],1:dims[2,j,i],j] <- object$results$estimation$paramLatentAuto[[j,i]]
			}
		}

		### Calculate the average coefficients
		paramLatentAutoMean <- vector("list",length=nAuto)
		for(i in 1:nAuto){
			paramLatentAutoMean[[i]]<-apply(paramLatentAuto[[i]], 1:2, mean)
		}

		### names each part of the object
		for(i in 1:nAuto){
			rownames(paramLatentAutoMean[[i]]) <- rownames(object$results$estimation$paramLatentAuto[[1,i]])
			colnames(paramLatentAutoMean[[i]]) <- colnames(object$results$estimation$paramLatentAuto[[which.max(dims[2,,i]),i]])
		}

		paramLatentAutoMean<-matrix(paramLatentAutoMean,ncol=nAuto)
		colnames(paramLatentAutoMean)<-colnames(object$results$estimation$paramLatentAuto)

		res$paramLatentAuto <- paramLatentAutoMean
	}

	### latentAuto
	if(any(names(res)=="latentAuto")){
		### Number of random effect
		nAuto <- ncol(object$results$estimation$latentAuto)

		### Number of latentAuto variables for each random effect
		dims <- array(NA, dim=c(2,niter, nAuto))

		for(i in 1:nAuto){
			dims[,,i] <- sapply(object$results$estimation$latentAuto[,i],dim)
		}

		### Store paramLatent for each random effect
		latentAuto <- vector("list",length=nAuto)
		for(i in 1:nAuto){
			latentAuto[[i]] <- array(0,dim=c(max(dims[1,,i]),max(dims[2,,i]),niter))
		}

		for(i in 1:nAuto){
			for(j in 1:niter){
				latentAuto[[i]][1:dims[1,j,i],1:dims[2,j,i],j] <- object$results$estimation$latentAuto[[j,i]]
			}
		}

		### Calculate the average coefficients
		latentAutoMean <- vector("list",length=nAuto)
		for(i in 1:nAuto){
			latentAutoMean[[i]]<-apply(latentAuto[[i]], 1:2, mean)
		}

		### names each part of the object
		for(i in 1:nAuto){
			rownames(latentAutoMean[[i]]) <- rownames(object$results$estimation$latentAuto[[1,i]])
			colnames(latentAutoMean[[i]]) <- colnames(object$results$estimation$latentAuto[[which.max(dims[2,,i]),i]])
		}

		latentAutoMean<-matrix(latentAutoMean,ncol=nAuto)
		colnames(latentAutoMean)<-colnames(object$results$estimation$latentAuto)

		res$latentAuto <- latentAutoMean
	}

	### paramAuto
	if(any(names(res)=="paramAuto")){
		### Number of random effect
		nAuto <- ncol(object$results$estimation$paramAuto)

		### Number of latentAuto variables for each random effect
		dims <- array(NA, dim=c(2,niter, nAuto))

		for(i in 1:nAuto){
			dims[,,i] <- sapply(object$results$estimation$paramAuto[,i],dim)
		}

		### Store paramLatent for each random effect
		paramAuto <- vector("list",length=nAuto)
		for(i in 1:nAuto){
			paramAuto[[i]] <- array(0,dim=c(max(dims[1,,i]),max(dims[2,,i]),niter))
		}

		for(i in 1:nAuto){
			for(j in 1:niter){
				paramAuto[[i]][1:dims[1,j,i],1:dims[2,j,i],j] <- object$results$estimation$paramAuto[[j,i]]
			}
		}

		### Calculate the average coefficients
		paramAutoMean <- vector("list",length=nAuto)
		for(i in 1:nAuto){
			paramAutoMean[[i]]<-apply(paramAuto[[i]], 1:2, mean)
		}

		### names each part of the object
		for(i in 1:nAuto){
			rownames(paramAutoMean[[i]]) <- rownames(object$results$estimation$paramAuto[[1,i]])
			colnames(paramAutoMean[[i]]) <- colnames(object$results$estimation$paramAuto[[which.max(dims[2,,i]),i]])
		}

		paramAutoMean<-matrix(paramAutoMean,ncol=nAuto)
		colnames(paramAutoMean)<-colnames(object$results$estimation$paramAuto)

		res$paramAuto <- paramAutoMean
	}
	return(res)
}
