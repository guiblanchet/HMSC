#' @title Full joint posterior distribution
#'
#' @description Calculate the full joint posterior distribution of the model for each iterations
#'
#'
#' @param hmsc An object of the class \code{hmsc}
#'
#' @return
#'
#' A three dimensional array with as first dimension the number of samples in the community matrix \code{Y}, the second dimension the number of species in the community matrix \code{Y} and as third dimension the number of iterations after the burning phase when constructing the model.
#'
#' @author F. Guillaume Blanchet
#'
#' @examples
#'
#' #================
#' ### Generate data
#' #================
#' desc <- cbind(scale(1:50), scale(1:50)^2)
#' nspecies <- 20
#' commDesc <- communitySimul(X = desc, nsp = nspecies)
#'
#' #=============
#' ### Formatting
#' #=============
#' ### Format data
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE, interceptTr = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 2000, nburn = 1000, thin = 1, verbose = 100)
#'
#' #==============================================
#' ### Calculate full joint posterior distribution
#' #==============================================
#' jpost <- jposterior(modelDesc)
#'
#' @keywords univar, multivariate, regression
#' @export
jposterior<-function(hmsc){

	### Basic objects
	nsite<-nrow(hmsc$data$Y)
	nsp<-ncol(hmsc$data$Y)
	nAuto<-length(hmsc$data$Auto)
	nRandom<-ncol(hmsc$data$Random)

	### Number of iterations
	if(!is.null(nAuto)){
		niter<-length(hmsc$results$estimation$paramLatentAuto)
	}
	if(!is.null(nRandom)){
		niter<-length(hmsc$results$estimation$paramLatent)
	}
	if(!is.null(hmsc$results$estimation$paramX)){
		niter<-dim(hmsc$results$estimation$paramX)[3]
	}

	### Results object
	res<-array(0,dim=c(nsite,nsp,niter))

	### Names each dimensions
	dimnames(res)[[1]]<-rownames(hmsc$data$Y)
	dimnames(res)[[2]]<-colnames(hmsc$data$Y)

	if(!is.null(nAuto)){
		dimnames(res)[[3]]<-rownames(hmsc$results$estimation$paramLatentAuto)
	}
	if(!is.null(nRandom)){
		dimnames(res)[[3]]<-rownames(hmsc$results$estimation$paramLatent)
	}
	if(!is.null(hmsc$results$estimation$paramX)){
		dimnames(res)[[3]]<-dimnames(hmsc$results$estimation$paramX)[[3]]
	}

	### Fill the results object
	if(is.null(hmsc$data$Random)){
		if(is.null(hmsc$data$Auto)){
			if(is.null(hmsc$data$X)){
				stop("This object is essentially empty")
			}else{
				### Only X
				for(i in 1:niter){
					res[,,i]<-tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
				}
			}
		}else{
			if(is.null(hmsc$data$X)){
				### Only Auto
				for(i in 1:niter){
					for(j in 1:nAuto){
						AutoModel<-tcrossprod(hmsc$results$estimation$latentAuto[[i,j]],hmsc$results$estimation$paramLatentAuto[[i,j]])
						res[,,i]<-res[,,i]+AutoModel[hmsc$data$Auto[[j]][,1],]
					}
				}
			}else{
				### X and Auto
				for(i in 1:niter){
					res[,,i]<-tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
					for(j in 1:nAuto){
						AutoModel<-tcrossprod(hmsc$results$estimation$latentAuto[[i,j]],hmsc$results$estimation$paramLatentAuto[[i,j]])
						res[,,i]<-res[,,i]+AutoModel[hmsc$data$Auto[[j]][,1],]
					}
				}
			}
		}
	}else{
		if(is.null(hmsc$data$Auto)){
			if(is.null(hmsc$data$X)){
				### Only Random
				for(i in 1:niter){
					for(j in 1:nRandom){
						RandomModel<-tcrossprod(hmsc$results$estimation$latent[[i,j]],hmsc$results$estimation$paramLatent[[i,j]])
						res[,,i]<-res[,,i]+RandomModel[hmsc$data$Random[,j],]
					}
				}
			}else{
				### X and Random
				for(i in 1:niter){
					res[,,i]<-tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
					for(j in 1:nRandom){
						RandomModel<-tcrossprod(hmsc$results$estimation$latent[[i,j]],hmsc$results$estimation$paramLatent[[i,j]])
						res[,,i]<-res[,,i]+RandomModel[hmsc$data$Random[,j],]
					}
				}
			}
		}else{
			if(is.null(hmsc$data$X)){
				### Auto and Random
				for(i in 1:niter){
					for(j in 1:nAuto){
						AutoModel<-tcrossprod(hmsc$results$estimation$latentAuto[[i,j]],hmsc$results$estimation$paramLatentAuto[[i,j]])
						res[,,i]<-res[,,i]+AutoModel[hmsc$data$Auto[[j]][,1],]
					}
					for(j in 1:nRandom){
						RandomModel<-tcrossprod(hmsc$results$estimation$latent[[i,j]],hmsc$results$estimation$paramLatent[[i,j]])
						res[,,i]<-res[,,i]+RandomModel[hmsc$data$Random[,j],]
					}
				}
			}else{
				### X, Auto and Random
				for(i in 1:niter){
					res[,,i]<-tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
					for(j in 1:nAuto){
						AutoModel<-tcrossprod(hmsc$results$estimation$latentAuto[[i,j]],hmsc$results$estimation$paramLatentAuto[[i,j]])
						res[,,i]<-res[,,i]+AutoModel[hmsc$data$Auto[[j]][,1],]
					}
					for(j in 1:nRandom){
						RandomModel<-tcrossprod(hmsc$results$estimation$latent[[i,j]],hmsc$results$estimation$paramLatent[[i,j]])
						res[,,i]<-res[,,i]+RandomModel[hmsc$data$Random[,j],]
					}
				}
			}
		}
	}

	return(res)
}
