#' @title Model prediction for HMSC models
#'
#' @description Obtains predictions for a model constructed using \code{hmsc}
#'
#' @param object An object of the class \code{hmsc}
#' @param newdata An optional object of class \code{HMSCdata} in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param type A character vector defining whether the prediction should be return on the scale of the response variable ("response", default) or on the scale of the linear predictors ("link"). See details.
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
#' #===================================
#' ### Calculate conditional prediction
#' #===================================
#' ### Conditional only on sp1
#' predCondSp1 <- predict(modelDesc, conditional = "sp1", nsample = 10)
#'
#' ### Conditional on the whole community but with an interest on sp1
#' formdata$Y[,1] <- NA
#' predCondAllSp1 <- predict(modelDesc, conditional = colnames(formdata$Y), nsample = 10)[,1]
#'
#' @keywords univar, multivariate, regression
#' @export
predict.hmsc<-function(object, newdata, type = "response", conditional = NULL, nsample, ...){

	type <- match.arg(type)

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
	if(is.null(conditional)){
		if(missing(newdata) || is.null(newdata)){
			if(is.null(data$Random)){
				if(is.null(data$Auto)){
					if(is.null(data$X)){
						stop("This object is essentially empty")
					}else{
						### Only X
						for(i in 1:niter){
							res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						}
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
						for(i in 1:niter){
							res[,,i]<-tcrossprod(data$X,object$results$estimation$paramX[,,i])
						}
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
	}

	### Fill the results object
	if(!is.null(conditional)){
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

		### Construct residVar for the different types of models
		if(any(class(object)=="probit") | any(class(object)=="poisson")){
			residVar <- matrix(1,nrow=niter,ncol=nsp)
		}

		if(any(class(object)=="gaussian")){
			residVar <- as.matrix(object$results$estimation$varNormal)[,spSep]
		}

		if(any(class(object)=="overPoisson")){
			residVar <- as.matrix(object$results$estimation$varPoisson)[,spSep]
		}

		if(is.null(data$Random)){
			if(is.null(data$Auto)){
				if(is.null(data$X)){
					stop("This object is essentially empty")
				}else{
					#=========
					### Only X
					#=========
					### Isolate parameters for selected species
					paramX <- array(object$results$estimation$paramX[spSel,,],dim=c(nsp,dim(object$results$estimation$paramX)[2],niter))

					res <- sampleCondPredX(Y,
						 										 data$X,
																 paramX,
																 residVar,
																 nsite,
																 nsp,
																 niter,
																 nsample,
																 class(object)[2])
				}
			}else{
				if(is.null(data$X)){
					#============
					### Only Auto
					#============
					if(!(missing(newdata) || is.null(newdata))){
						latentAuto<-vector("list",length=niter*nAuto)
						for(i in 1:niter){
							for(j in 1:nAuto){
								### Construct latentAuto variables
								nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
								latentAuto[[i,j]]<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
								rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])
							}
						}
					}else{
						latentAuto <- object$results$estimation$latentAuto
					}

					### Build Auto and RandomAuto
					RandomAuto <- data.frame(lapply(data$Auto,function(x) x[,1]))
					nAutoLev <- sapply(RandomAuto, nlevels)
					RandomAuto<-sapply(RandomAuto,as.numeric)-1

					AutoCoord<-vector("list",length=nAuto)

					for(i in 1:nAuto){
						nAutoCoord<-ncol(data$Auto[[i]])-1
						AutoCoordMean<-matrix(NA,nrow=nAutoLev[i],ncol=nAutoCoord)

						for(j in 1:nAutoCoord){
							AutoCoordMean[,j]<-tapply(data$Auto[[i]][,j+1],data$Auto[[i]][,1],mean)
						}

						AutoCoord[[i]]<-AutoCoordMean
					}

					### Flat prior definition
					priorParamAutoDist <- flatPriorAuto(data, family = class(object)[2])$paramAutoDist
					npriorParamAuto <- nrow(priorParamAutoDist)

					### Isolate parameters for selected species
					paramLatentAuto<-vector("list",length=niter*nAuto)
					dim(paramLatentAuto)<-c(niter,nAuto)
					for(i in 1:niter){
						for(j in 1:nAuto){
							paramLatentAuto[[i,j]] <- matrix(object$results$estimation$paramLatentAuto[[i,j]][spSel,],nrow=nsp)
						}
					}

					res <- sampleCondPredAuto(Y,
																		AutoCoord,
																		RandomAuto,
																		latentAuto,
																		paramLatentAuto,
																		object$results$estimation$paramAuto,
																		residVar,
																		priorParamAutoDist,
																		nsite,
																		nsp,
																		nAuto,
																		nAutoLev,
																		npriorParamAuto,
																		niter,
																		nsample,
																		class(object)[2])
				}else{
					#=============
					### X and Auto
					#=============
					if(!(missing(newdata) || is.null(newdata))){
						latentAuto<-vector("list",length=niter*nAuto)
						for(i in 1:niter){
							for(j in 1:nAuto){
								### Construct latentAuto variables
								nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
								latentAuto[[i,j]]<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
								rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])
							}
						}
					}else{
						latentAuto <- object$results$estimation$latentAuto
					}

					### Build Auto and RandomAuto
					Auto <- lapply(data$Auto,function(x) as.matrix(x[,-1]))
					RandomAuto <- data.frame(lapply(data$Auto,function(x) x[,1]))
					nAutoLev <- sapply(RandomAuto, nlevels)
					RandomAuto<-sapply(RandomAuto,as.numeric)-1

					### Flat prior definition
					priorParamAutoDist <- flatPriorAuto(data, family = class(object)[2])$paramAutoDist
					npriorParamAuto <- nrow(priorParamAutoDist)

					### Isolate parameters for selected species
					paramX <- array(object$results$estimation$paramX[spSel,,],dim=c(nsp,dim(object$results$estimation$paramX)[2],niter))

					paramLatentAuto<-vector("list",length=niter*nAuto)
					dim(paramLatentAuto)<-c(niter,nAuto)
					for(i in 1:niter){
						for(j in 1:nAuto){
							paramLatentAuto[[i,j]] <- matrix(object$results$estimation$paramLatentAuto[[i,j]][spSel,],nrow=nsp)
						}
					}

					res <- sampleCondPredXAuto(Y,
																		data$X,
																		Auto,
																		RandomAuto,
																		paramX,
																		latentAuto,
																		paramLatentAuto,
																		object$results$estimation$paramAuto,
																		residVar,
																		priorParamAutoDist,
																		nsite,
																		nsp,
																		nAuto,
																		nAutoLev,
																		npriorParamAuto,
																		niter,
																		nsample,
																		class(object)[2])
				}
			}
		}else{
			if(is.null(data$Auto)){
				if(is.null(data$X)){
					#==============
					### Only Random
					#==============
					if(!(missing(newdata) || is.null(newdata))){
						latent<-vector("list",length=niter*nRandom)
						dim(latent)<-c(niter,nRandom)
						for(i in 1:niter){
							for(j in 1:nRandom){
								### Construct latent variables
								nlatent<-ncol(object$results$estimation$latent[[i,j]])
								latent[[i,j]] <- rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
								rownames(latent[[i,j]])<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])

								### Isolate parameters for selected species
								paramLatent<-vector("list",length=niter*nRandom)
								dim(paramLatent)<-c(niter,nRandom)
								paramLatent[[i,j]] <- object$results$estimation$paramLatent[[i,j]][spSel,]
							}
						}
					}else{
						latent <- object$results$estimation$latent
					}

					### Isolate parameters for selected species
					paramLatent<-vector("list",length=niter*nRandom)
					dim(paramLatent)<-c(niter,nRandom)
					for(i in 1:niter){
						for(j in 1:nRandom){
							paramLatent[[i,j]] <- matrix(object$results$estimation$paramLatent[[i,j]][spSel,],nrow=nsp)
						}
					}

					### Organize Random
					Random<-sapply(data$Random,as.numeric)-1
					nRandomLev <- sapply(data$Random, nlevels)

					res <- sampleCondPredLatent(Y,
																			Random,
																			latent,
																			paramLatent,
																			residVar,
																			nsite,
																			nsp,
																			nRandom,
																			nRandomLev,
																			niter,
																			nsample,
																			class(object)[2])
				}else{
					#===============
					### X and Random
					#===============
					if(!(missing(newdata) || is.null(newdata))){
						latent<-vector("list",length=niter*nRandom)
						dim(latent)<-c(niter,nRandom)
						for(i in 1:niter){
							for(j in 1:nRandom){
								### Construct latent variables
								nlatent<-ncol(object$results$estimation$latent[[i,j]])
								latent[[i,j]] <- rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
								rownames(latent[[i,j]])<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])
						 	}
						}
					}else{
						latent <- object$results$estimation$latent
					}

					### Organize Random
					Random<-sapply(data$Random,as.numeric)-1
					nRandomLev <- sapply(data$Random, nlevels)

					### Isolate parameters for selected species
					paramX <- array(object$results$estimation$paramX[spSel,,],dim=c(nsp,dim(object$results$estimation$paramX)[2],niter))

					paramLatent<-vector("list",length=niter*nRandom)
					dim(paramLatent)<-c(niter,nRandom)
					for(i in 1:niter){
						for(j in 1:nRandom){
							paramLatent[[i,j]] <- matrix(object$results$estimation$paramLatent[[i,j]][spSel,],nrow=nsp)
						}
					}

					res <- sampleCondPredXLatent(Y,
						 													 data$X,
																			 Random,
																			 paramX,
																			 latent,
																			 paramLatent,
																			 residVar,
																			 nsite,
																			 nsp,
																			 nRandom,
																			 nRandomLev,
																			 niter,
																			 nsample,
																			 class(object)[2])
				}
			}else{
				if(is.null(data$X)){
					#==================
					### Auto and Random
					#==================
					if(!(missing(newdata) || is.null(newdata))){
						latent<-vector("list",length=niter*nRandom)
						latentAuto<-vector("list",length=niter*nAuto)
						dim(latent)<-c(niter,nRandom)
						for(i in 1:niter){
							for(j in 1:nRandom){
								### Construct latent variables
								nlatent<-ncol(object$results$estimation$latent[[i,j]])
								latent[[i,j]] <- rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
								rownames(latent[[i,j]])<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])
						 	}

							for(j in 1:nAuto){
								### Construct latentAuto variables
								nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
								latentAuto[[i,j]]<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
								rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])
							}
						}
					}else{
						latent <- object$results$estimation$latent
						latentAuto <- object$results$estimation$latentAuto
					}

					### Organize Random
					Random<-sapply(data$Random,as.numeric)-1
					nRandomLev <- sapply(data$Random, nlevels)

					### Build Auto and RandomAuto
					Auto <- lapply(data$Auto,function(x) as.matrix(x[,-1]))
					RandomAuto <- data.frame(lapply(data$Auto,function(x) x[,1]))
					nAutoLev <- sapply(RandomAuto, nlevels)
					RandomAuto<-sapply(RandomAuto,as.numeric)-1

					### Flat prior definition
					priorParamAutoDist <- flatPriorAuto(data, family = class(object)[2])$paramAutoDist
					npriorParamAuto <- nrow(priorParamAutoDist)

					### Isolate parameters for selected species
					paramLatent<-vector("list",length=niter*nRandom)
					dim(paramLatent)<-c(niter,nRandom)
					for(i in 1:niter){
						for(j in 1:nRandom){
							paramLatent[[i,j]] <- matrix(object$results$estimation$paramLatent[[i,j]][spSel,],nrow=nsp)
						}
					}

					paramLatentAuto<-vector("list",length=niter*nAuto)
					dim(paramLatentAuto)<-c(niter,nAuto)
					for(i in 1:niter){
						for(j in 1:nAuto){
							paramLatentAuto[[i,j]] <- matrix(object$results$estimation$paramLatentAuto[[i,j]][spSel,],nrow=nsp)
						}
					}

					res <- sampleCondPredLatentAuto(Y,
																					Auto,
																					data$Random,
																					RandomAuto,
																					latent,
																					paramLatent,
																					latentAuto,
																					paramLatentAuto,
																					object$results$estimation$paramAuto,
																					residVar,
																					priorParamAutoDist,
																					nsite,
																					nsp,
																					nRandom,
																					nRandomLev,
																					nAuto,
																					nAutoLev,
																					npriorParamAuto,
																					niter,
																					nsample,
																					class(object)[2])
				}else{
					#=====================
					### X, Auto and Random
					#=====================
					if(!(missing(newdata) || is.null(newdata))){
						latent<-vector("list",length=niter*nRandom)
						latentAuto<-vector("list",length=niter*nAuto)
						dim(latent)<-c(niter,nRandom)
						for(i in 1:niter){
							for(j in 1:nRandom){
								### Construct latent variables
								nlatent<-ncol(object$results$estimation$latent[[i,j]])
								latent[[i,j]] <- rbind(object$results$estimation$latent[[i,j]][matchRandom[[j]],],matrix(rnorm(length(noMatchRandom[[j]])*nlatent),ncol=nlatent))
								rownames(latent[[i,j]])<-c(objectRandomNlev[[j]][matchRandom[[j]]],noMatchRandom[[j]])
						 	}

							for(j in 1:nAuto){
								### Construct latentAuto variables
								nlatentAuto<-ncol(object$results$estimation$latentAuto[[i,j]])
								latentAuto[[i,j]]<-rbind(object$results$estimation$latentAuto[[i,j]][matchAuto[[j]],],matrix(rnorm(length(noMatchAuto[[j]])*nlatentAuto),ncol=nlatentAuto))
								rownames(latentAuto)<-c(objectAutoNlev[[j]][matchAuto[[j]]],noMatchAuto[[j]])
							}
						}
					}else{
						latent <- object$results$estimation$latent
						latentAuto <- object$results$estimation$latentAuto
					}

					### Organize Random
					Random<-sapply(data$Random,as.numeric)-1
					nRandomLev <- sapply(data$Random, nlevels)

					### Build Auto and RandomAuto
					Auto <- lapply(data$Auto,function(x) as.matrix(x[,-1]))
					RandomAuto <- data.frame(lapply(data$Auto,function(x) x[,1]))
					nAutoLev <- sapply(RandomAuto, nlevels)
					RandomAuto<-sapply(RandomAuto,as.numeric)-1

					### Flat prior definition
					priorParamAutoDist <- flatPriorAuto(data, family = class(object)[2])$paramAutoDist
					npriorParamAuto <- nrow(priorParamAutoDist)

					### Isolate parameters for selected species
					paramX <- array(object$results$estimation$paramX[spSel,,],dim=c(nsp,dim(object$results$estimation$paramX)[2],niter))

					paramLatent<-vector("list",length=niter*nRandom)
					dim(paramLatent)<-c(niter,nRandom)
					for(i in 1:niter){
						for(j in 1:nRandom){
							paramLatent[[i,j]] <- matrix(object$results$estimation$paramLatent[[i,j]][spSel,],nrow=nsp)
						}
					}

					paramLatentAuto<-vector("list",length=niter*nAuto)
					dim(paramLatentAuto)<-c(niter,nAuto)
					for(i in 1:niter){
						for(j in 1:nAuto){
							paramLatentAuto[[i,j]] <- matrix(object$results$estimation$paramLatentAuto[[i,j]][spSel,],nrow=nsp)
						}
					}

					res <- sampleCondPredXLatentAuto(Y,
																					 data$X,
																					 Auto,
																					 data$Random,
																					 RandomAuto,
																					 paramX,
																					 latent,
																					 paramLatent,
																					 latentAuto,
																					 paramLatentAuto,
																					 object$results$estimation$paramAuto,
																					 residVar,
																					 priorParamAutoDist,
																					 nsite,
																					 nsp,
																					 nRandom,
																					 nRandomLev,
																					 nAuto,
																					 nAutoLev,
																					 npriorParamAuto,
																					 niter,
																					 nsample,
																					 class(object)[2])
				}
			}
		}

		### Average the result matrix
		res <- apply(res,1:2, mean)
	}else{
		Y <- data$Y
	}

	if(is.null(conditional)){
		### Apply inverse link function
		if(type=="response"){
			if(any(class(object)=="probit")){
				result<-pnorm(apply(res,1:2, mean))
			}

			if(any(class(object)=="gaussian")){
				result<-apply(res,1:2, mean)
			}

			if(any(class(object)=="poisson" | any(class(object)=="overPoisson"))){
				result<-exp(apply(res,1:2, mean))
			}
		}

		if(type=="link"){
				result<-apply(res,1:2, mean)
		}
		colnames(result)<-colnames(Y)
	}else{
		### Apply inverse link function
		if(type=="response"){
			if(any(class(object)=="probit")){
				result<-pnorm(res)
			}

			if(any(class(object)=="gaussian")){
				result<-res
			}

			if(any(class(object)=="poisson" | any(class(object)=="overPoisson"))){
				result<-exp(res)
			}
		}

		if(type=="link"){
			result<-res
		}
		colnames(result)<-colnames(Y)
	}

	### Return model
	return(result)
}
