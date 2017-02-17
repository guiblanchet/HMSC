#' @title Model prediction for HMSC models
#'
#' @description Obtains predictions for a model constructed using \code{hmsc}
#'
#' @param object An object of the class \code{hmsc}
#' @param newdata An optional of class \code{HMSCdata} in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param \dots Additional arguments affecting the predictions produced.
#' 
#' @return
#'
#' A matrix with the predicted values for each species at the sites in \code{newdata} or a matrix of the same dimension as the fitted values. 
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
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE, 
#' 					     interceptTr = FALSE)
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
predict.hmsc<-function(object,newdata,...){
	
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
	
	### Apply inverse link function
	if(any(class(object)=="probit")){
#		result<-list(iter=pnorm(res),mean=pnorm(apply(res,1:2, mean)))
		result<-pnorm(apply(res,1:2, mean))
	}
	
	if(any(class(object)=="gaussian")){
#		result<-list(iter=res,mean=apply(res,1:2, mean))
		result<-apply(res,1:2, mean)
	}
	
#	colnames(result$mean)<-colnames(object$data$Y)
	colnames(result)<-colnames(object$data$Y)
	
	### Return model
	return(result)
}
