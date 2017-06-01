#' @title Model prediction for HMSC models
#'
#' @description Obtains predictions for a model constructed using \code{hmsc}
#'
#' @param object An object of the class \code{hmsc}
#' @param newdata An optional object of class \code{HMSCdata} in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param conditional A character vector defining the names of the species used for the conditional prediction. If omitted, unconditional predictions are carried out. Default is NULL.
#' @param nsample A numerical value defining the number of samples to carry out when calculating conditional probability. If the \code{conditional} is NULL, this argument will not be considered.
#' @param \dots Additional arguments affecting the predictions produced.
#'
#' @details
#'
#' The function is designed for conditional predictions to be calculated so species are the variables that can be used to condition the prediction.
#'
#' @return
#'
#' A matrix with the predicted values for each species at the sites in \code{newdata} or a matrix of the same dimension as the fitted values.
#'
#' If conditional predictions are performed, the result matrix is only given for the species considered.
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
#' #=======================
#' ### Calculate prediction
#' #=======================
#' predModel <- predict(modelDesc)
#'
#' @keywords univar, multivariate, regression
#' @export
predict.hmsc<-function(object, newdata, conditional=NULL, nsample, ...){

	### Data to use for prediction
	if(missing(newdata) || is.null(newdata)){
		data<-object$data
	}else{
		### A few checks
		if(class(newdata)!="HMSCdata"){
			stop("'newdata' needs to be of class 'HMSCdata'")
		}else{
			### Check if names match
			namesNewdata<-names(newdata)
			namesObject<-names(object$data)

			if(any(namesNewdata=="Y")){
				if(length(namesNewdata)!=length(namesObject)){
					stop("'newdata' has a different number of datasets than 'object$data' and cannot be used for prediction")
				}
			}else{
				namesObject<-names(object$data)[-1] # Remove Y
				if(length(namesNewdata)!=length(namesObject)){
					stop("'newdata' has a different number of datasets than 'object$data' and cannot be used for prediction")
				}
			}

			if(!all(namesNewdata==namesObject)){
				stop("the data in 'newdata' has different datasets than the ones present in 'object$data' and cannot be used for prediction")
			}else{
				### Check to make sure that all the levels in newdata$Auto are present in object$data$Auto
				if(!is.null(object$data$Auto)){
					if(length(object$data$Auto)!=length(newdata$Auto)){
						stop("'newdata' does not have the same number of autocorrelated factors as 'object$data'")
					}else{
						objectAutoNlev<-lapply(object$data$Auto,function(x) levels(x[,1]))
						newdataAutoNlev<-lapply(newdata$Auto,function(x) levels(x[,1]))

						### Find the levels that are in both objectAutoNlev and newdataAutoNlev
						matchAuto<-vector("list",length=length(objectAutoNlev)) # List of values (pointers)
						noMatchAuto<-vector("list",length=length(objectAutoNlev)) # List of names

						for(i in 1:length(objectAutoNlev)){
							matchAuto[[i]]<-which(objectAutoNlev[[i]] %in% newdataAutoNlev[[i]])

							if(length(matchAuto[[i]])>0){
								noMatchAuto[[i]]<-newdataAutoNlev[[i]][-matchAuto[[i]]]
							}else{
								noMatchAuto[[i]]<-newdataAutoNlev[[i]]
							}
						}
					}
				}

				### Check to make sure that all the levels in newdata$Random are present in object$data$Random
				if(!is.null(object$data$Random)){
					if(ncol(object$data$Random)!=ncol(newdata$Random)){
						stop("'newdata' does not have the same number of random factors as 'object$data'")
					}else{
						objectRandomNlev<-lapply(object$data$Random, levels)
						newdataRandomNlev<-lapply(newdata$Random,levels)

						### Find the levels that are in both objectRandomNlev and newdataRandomNlev
						matchRandom<-vector("list",length=length(objectRandomNlev)) # List of values (pointers)
						noMatchRandom<-vector("list",length=length(objectRandomNlev)) # List of names

						for(i in 1:length(objectRandomNlev)){
							matchRandom[[i]]<-which(objectRandomNlev[[i]] %in% newdataRandomNlev[[i]])

							if(length(matchRandom[[i]])>0){
								noMatchRandom[[i]]<-newdataRandomNlev[[i]][-matchRandom[[i]]]
							}else{
								noMatchRandom[[i]]<-newdataRandomNlev[[i]]
							}
						}
					}
				}

				### Check to make sure that all the variables in newdata$X  are present in object$data$X
				if(!is.null(object$data$X)){
					if(ncol(object$data$X)!=ncol(newdata$X)){
						stop("The number of X variables is different in 'newdata' and 'object$data'")
					}
					if(!all(colnames(object$data$X)==colnames(newdata$X))){
						stop("Some X variables in 'newdata' are not in 'object$data'")
					}
				}
			}
		}

		### New data used for prediction
		data<-newdata
	}

	### Basic objects
	nAuto<-length(data$Auto)
	nRandom<-ncol(data$Random)
	nsp<-ncol(object$data$Y)

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

	### Number of sites
	if(!is.null(nAuto)){
		nsite<-nrow(data$Auto[[1]])
	}
	if(!is.null(nRandom)){
		nsite<-nrow(data$Random)
	}
	if(!is.null(object$results$estimation$paramX)){
		nsite<-nrow(data$X)
	}

	### Results object
	res<-array(0,dim=c(nsite,nsp,niter))

	### Names each dimensions
	dimnames(res)[[1]]<-rownames(data$Y)
	dimnames(res)[[2]]<-colnames(data$Y)

	if(!is.null(nAuto)){
		dimnames(res)[[3]]<-rownames(object$results$estimation$paramLatentAuto)
	}
	if(!is.null(nRandom)){
		dimnames(res)[[3]]<-rownames(object$results$estimation$paramLatent)
	}
	if(!is.null(object$results$estimation$paramX)){
		dimnames(res)[[3]]<-dimnames(object$results$estimation$paramX)[[3]]
	}

	### Fill the results object
	if(missing(newdata) || is.null(newdata)){
		if(is.null(data$Random)){
			if(is.null(data$Auto)){
				if(is.null(data$X)){
					stop("This object is essentially empty")
				}else{
					### Only X
					res<-predLinkX(data$X,object$results$estimation$paramX,nsite,nsp,niter)
				}
			}else{
				if(is.null(data$X)){
					### Only Auto
					for(i in 1:niter){
						for(j in 1:nAuto){
							AutoModel<-tcrossprod(object$results$estimation$latentAuto[[i,j]],object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
					}
				}else{
					### X and Auto
					for(i in 1:niter){
						res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						for(j in 1:nAuto){
							AutoModel<-tcrossprod(object$results$estimation$latentAuto[[i,j]],object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
					}
				}
			}
		}else{
			if(is.null(data$Auto)){
				if(is.null(data$X)){
					### Only Random
					for(i in 1:niter){
						for(j in 1:nRandom){
							RandomModel<-tcrossprod(object$results$estimation$latent[[i,j]],object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}else{
					### X and Random
					for(i in 1:niter){
						res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						for(j in 1:nRandom){
							RandomModel<-tcrossprod(object$results$estimation$latent[[i,j]],object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}
			}else{
				if(is.null(data$X)){
					### Auto and Random
					for(i in 1:niter){
						for(j in 1:nAuto){
							AutoModel<-tcrossprod(object$results$estimation$latentAuto[[i,j]],object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
						for(j in 1:nRandom){
							RandomModel<-tcrossprod(object$results$estimation$latent[[i,j]],object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}else{
					### X, Auto and Random
					for(i in 1:niter){
						res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						for(j in 1:nAuto){
							AutoModel<-tcrossprod(object$results$estimation$latentAuto[[i,j]],object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
						for(j in 1:nRandom){
							RandomModel<-tcrossprod(object$results$estimation$latent[[i,j]],object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}
			}
		}
	}else{
		if(is.null(data$Random)){
			if(is.null(data$Auto)){
				if(is.null(data$X)){
					stop("This object is essentially empty")
				}else{
					### Only X
					res<-predLinkX(data$X,object$results$estimation$paramX,nsite,nsp,niter)
				}
			}else{
				if(is.null(data$X)){
					### Only Auto
					for(i in 1:niter){
						for(j in 1:nAuto){
							### Construct latentAuto variables
							nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
							latentAuto<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
							rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])

							AutoModel<-tcrossprod(latentAuto,object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
					}
				}else{
					### X and Auto
					for(i in 1:niter){
						res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						for(j in 1:nAuto){
							### Construct latentAuto variables
							nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
							latentAuto<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
							rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])

							AutoModel<-tcrossprod(latentAuto,object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
					}
				}
			}
		}else{
			if(is.null(data$Auto)){
				if(is.null(data$X)){
					### Only Random
					for(i in 1:niter){
						for(j in 1:nRandom){
							### Construct latent variables
							nlatent<-ncol(object$results$estimation$latent[[i,j]])
							latent<-rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
							rownames(latent)<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])

							RandomModel<-tcrossprod(latent,object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}else{
					### X and Random
					for(i in 1:niter){
						res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						for(j in 1:nRandom){
							### Construct latent variables
							nlatent<-ncol(object$results$estimation$latent[[i,j]])
							latent<-rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
							rownames(latent)<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])

							RandomModel<-tcrossprod(latent,object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}
			}else{
				if(is.null(data$X)){
					### Auto and Random
					for(i in 1:niter){
						for(j in 1:nAuto){
							### Construct latentAuto variables
							nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
							latentAuto<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
							rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])

							AutoModel<-tcrossprod(latentAuto,object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
						for(j in 1:nRandom){
							### Construct latent variables
							nlatent<-ncol(object$results$estimation$latent[[i,j]])
							latent<-rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
							rownames(latent)<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])

							RandomModel<-tcrossprod(latent,object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}else{
					### X, Auto and Random
					for(i in 1:niter){
						res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						for(j in 1:nAuto){
							### Construct latentAuto variables
							nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
							latentAuto<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
							rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])

							AutoModel<-tcrossprod(latentAuto,object$results$estimation$paramLatentAuto[[i,j]])
							res[,,i]<-res[,,i]+AutoModel[data$Auto[[j]][,1],]
						}
						for(j in 1:nRandom){
							### Construct latent variables
							nlatent<-ncol(object$results$estimation$latent[[i,j]])
							latent<-rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
							rownames(latent)<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])

							RandomModel<-tcrossprod(latent,object$results$estimation$paramLatent[[i,j]])
							res[,,i]<-res[,,i]+RandomModel[data$Random[,j],]
						}
					}
				}
			}
		}
	}


	### Fill the results object
	if(!is.null(conditional)){
		### Average the result matrix
		res <- apply(res,1:2, mean)

		### A few checks
		if(!is.character(conditional)){
			stop("'conditional' needs to be a vector of characters")
		}

		### Find the species to consider as conditional
		condSp <- colnames(data$Y)%in%conditional
		if(sum(condSp)!=length(conditional)){
			stop("Some species defined in 'conditional' were not found in the species data matrix")
		}
		### Selected species
		spSel <- which(condSp)

		### Select the species in Y
		Y <- as.matrix(data$Y[,spSel])

		### Redefine the number of species
		nsp <- ncol(Y)

		### Extract the species to consider in the estimated model calculated above
		EstModel <- res[,spSel]

		### Check to make sure that EstModel always has 2 dimensions
		if(!is.matrix(EstModel)){
			EstModel <- as.matrix(EstModel)
		}

		### Construct residVar for the different types of models
		if(any(class(object)=="probit") | any(class(object)=="poisson")){
			residVar <- matrix(1,nrow=niter,ncol=nsp)
		}

		if(any(class(object)=="gaussian")){
			residVar <- as.matrix(object$results$estimation$varNormal[,spSel])
		}

		if(any(class(object)=="overPoisson")){
			residVar <- as.matrix(object$results$estimation$varPoisson[,spSel])
		}

		### Sample conditional prediction
		res <- sampleCondPred(Y, EstModel, residVar, nsite, nsp, nsample, family=class(object)[2])
	}else{
		Y <- data$Y
	}

	### Apply inverse link function
	if(any(class(object)=="probit")){
		result<-pnorm(apply(res,1:2, mean))
	}

	if(any(class(object)=="gaussian")){
		result<-apply(res,1:2, mean)
	}

	if(any(class(object)=="poisson" | any(class(object)=="overPoisson"))){
		result<-exp(apply(res,1:2, mean))
	}

	colnames(result)<-colnames(Y)

	### Return model
	return(result)
}
