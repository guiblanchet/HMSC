#' @title Variation partitioning 
#'
#' @description Partition the variation of models 
#'
#' @param hmsc An object of the class \code{hmsc}.
#' @param groupX A vector defining how the covariates (\code{X}) are grouped. This argument is ignored if the models does not have any covariates.
#' 
#' @return
#'
#' A list presenting how the explained variation in the data is partitioned between the different groups of the covariates (first part of thet list) and the random effects (second part of the list). Note that both autocorrelated (\code{Auto}) and non-autocorrelated (\code{Random}) are considered here. In addition, this function also calculate the amount of variation explained by traits.
#'
#' @author Gleb Tikhonov and modified by Guillaume Blanchet
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
#' 						   interceptTr = FALSE)
#' 
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 2000, nburn = 1000, thin = 1, verbose = 100)
#' 
#' #==================================
#' ### Calculate variance partitioning
#' #==================================
#' vp <- variPart(modelDesc, groupX = c("A", "B"))
#' vp
#'
#' #===================
#' ### Plot the results
#' #===================
#' barplot(t(vp), las = 2)
#' 
#' @keywords univar, multivariate, regression
#' @export
variPart<-function(hmsc,groupX){
	
	### Basic objects
	nsite<-nrow(hmsc$data$Y)
	nsp<-ncol(hmsc$data$Y)
	nX<-length(hmsc$data$X)
	nAuto<-length(hmsc$data$Auto)
	nRandom<-ncol(hmsc$data$Random)
	
	### Objects related to X
	if(!is.null(nX)){
		if(nX>0){
			nGroups<-length(unique(groupX))
			groupXLev<-unique(groupX)
			covX<-cov(hmsc$data$X)
			
			### Result objects
			variPartX <- vector(length=nsp,"numeric")
			variPartXSplit <- matrix(0, nrow=nsp, ncol=nGroups)
			rownames(variPartXSplit)<-colnames(hmsc$data$Y)
			colnames(variPartXSplit)<-groupXLev
			
			if(any(names(hmsc$data) == "Tr")){
				traitR2 <- 0
			}
		}
	}
	
	### Objects related to Auto
	if(!is.null(nAuto)){
		if(nAuto>0){
			variPartAuto <- matrix(0, nrow=nsp, ncol=nAuto)
			rownames(variPartAuto)<-colnames(hmsc$data$Y)
			colnames(variPartAuto)<-names(hmsc$data$Auto)
		}
	}
	
	### Objects related to Random
	if(!is.null(nRandom)){
		if(nRandom>0){
			variPartRandom <- matrix(0, nrow=nsp, ncol=nRandom)
			rownames(variPartRandom)<-colnames(hmsc$data$Y)
			colnames(variPartRandom)<-colnames(hmsc$data$Random)
		}
	}
	
	### Number of iterations
	if(!is.null(nAuto)){
		niter<-length(hmsc$results$estimation$paramLatentAuto)
	}
	if(!is.null(nRandom)){
		niter<-length(hmsc$results$estimation$paramLatent)
	}
	if(!is.null(nX)){
		niter<-dim(hmsc$results$estimation$paramX)[3]
	}
	
	### Fill the results object
	if(is.null(hmsc$data$Random)){
		if(is.null(hmsc$data$Auto)){
			if(is.null(hmsc$data$X)){
				stop("This object is essentially empty")
			}else{
				### Only X
				for(i in 1:niter){
					predXTotal <- rep(0, nsp)
					predXSplit <- matrix(0, nrow=nsp, ncol=nGroups)
					
					if(any(names(hmsc$data) == "Tr")){
						### Calculate R2 related to traits
						predX <- tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
						
						meansParamX <- t(hmsc$results$estimation$paramTr[,,i]%*%hmsc$data$Tr)
						predTr <- tcrossprod(hmsc$data$X,meansParamX)
						traitR2Base <- vector(length=2,"numeric")
						
						for(j in 1:nsite){
							covPred <- cov(cbind(predTr[j,], predX[j,]))
							traitR2Base[1] <- traitR2Base[1] + covPred[1,2]*covPred[2,1]
							traitR2Base[2] <- traitR2Base[2] + covPred[1,1]*covPred[2,2]
						}
						
						traitR2 <- traitR2+traitR2Base[1]/traitR2Base[2]
					}
					
					
					### Calculate R2 related to X
					for(j in 1:nsp){
						predXTotalSub <- hmsc$results$estimation$paramX[j,,i] %*% crossprod(covX,hmsc$results$estimation$paramX[j,,i])
						predXTotal[j] <- predXTotal[j] + predXTotalSub
						for(k in 1:nGroups){
							sel <- groupX==groupXLev[k]
							predXPart <- hmsc$results$estimation$paramX[j,sel,i] %*% crossprod(covX[sel,sel],hmsc$results$estimation$paramX[j,sel,i])
							predXSplit[j,k] <- predXSplit[j,k] + predXPart
						}
					}
					
					variPartX <- variPartX + rep(1, nsp)
					variPartXSplit <- variPartXSplit + predXSplit/replicate(nGroups, rowSums(predXSplit))
				}
				
				
				variPartX <- variPartX / niter
				variPartXSplit <- variPartXSplit / niter
				variPart<-replicate(nGroups,variPartX)*variPartXSplit;
				
				if(any(names(hmsc$data) == "Tr")){
					traitR2 <- traitR2 / niter
					res<-list(variPart=variPart,
							  traitR2=traitR2)
				}else{
					res<-variPart
				}
			}
		}else{
			if(is.null(hmsc$data$X)){
				### Only Auto
				for(i in 1:niter){
					PredAuto <- matrix(0, nrow=nsp, ncol=nAuto)
					for(j in 1:nAuto){
						paramLatentAuto<-hmsc$results$estimation$paramLatentAuto[[i,j]]
						nLatentAuto<-ncol(paramLatentAuto)
						
						for(k in 1:nLatentAuto){
							PredAuto[,j] <- PredAuto[,j] + paramLatentAuto[,k]*paramLatentAuto[,k]
						}
					}
				
					variPartAuto <- variPartAuto + PredAuto/replicate(nAuto, rowSums(PredAuto))
				}
				
				variPartAuto <- variPartAuto / niter
				res<-variPartAuto
			}else{
				### X and Auto
				for(i in 1:niter){
					predXTotal <- rep(0, nsp)
					predXSplit <- matrix(0, nrow=nsp, ncol=nGroups)
					PredAuto <- matrix(0, nrow=nsp, ncol=nAuto)
					
					if(any(names(hmsc$data) == "Tr")){
						### Calculate R2 related to traits
						predX <- tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
						
						meansParamX <- t(hmsc$results$estimation$paramTr[,,i]%*%hmsc$data$Tr)
						predTr <- tcrossprod(hmsc$data$X,meansParamX)
						traitR2Base <- vector(length=2,"numeric")
						
						for(j in 1:nsite){
							covPred <- cov(cbind(predTr[j,], predX[j,]))
							traitR2Base[1] <- traitR2Base[1] + covPred[1,2]*covPred[2,1]
							traitR2Base[2] <- traitR2Base[2] + covPred[1,1]*covPred[2,2]
						}
						
						traitR2 <- traitR2+traitR2Base[1]/traitR2Base[2]
					}
					
					
					### Calculate R2 related to X
					for(j in 1:nsp){
						predXTotalSub <- hmsc$results$estimation$paramX[j,,i] %*% crossprod(covX,hmsc$results$estimation$paramX[j,,i])
						predXTotal[j] <- predXTotal[j] + predXTotalSub
						for(k in 1:nGroups){
							sel <- groupX==groupXLev[k]
							predXPart <- hmsc$results$estimation$paramX[j,sel,i] %*% crossprod(covX[sel,sel],hmsc$results$estimation$paramX[j,sel,i])
							predXSplit[j,k] <- predXSplit[j,k] + predXPart
						}
					}
					
					### Calculate R2 related to auto
					for(j in 1:nAuto){
						paramLatentAuto<-hmsc$results$estimation$paramLatentAuto[[i,j]]
						nLatentAuto<-ncol(paramLatentAuto)
						
						for(k in 1:nLatentAuto){
							PredAuto[,j] <- PredAuto[,j] + paramLatentAuto[,k]*paramLatentAuto[,k]
						}
					}
				
					variTotal <- predXTotal + rowSums(PredAuto)
					variPartX <- variPartX + predXTotal/variTotal;
					variPartAuto <- variPartAuto + PredAuto/replicate(nAuto, variTotal)
					variPartXSplit <- variPartXSplit + predXSplit/replicate(nGroups, apply(predXSplit, 1, sum));
				}
				
				
				variPartX <- variPartX / niter
				variPartAuto <- variPartAuto / niter
				variPartXSplit <- variPartXSplit / niter
				variPart<-cbind(replicate(nGroups,variPartX)*variPartXSplit,variPartAuto)

				if(any(names(hmsc$data) == "Tr")){
					traitR2 <- traitR2 / niter
					res<-list(variPart=variPart,
							  traitR2=traitR2)
				}else{
					res<-variPart
				}
			}
		}
	}else{
		if(is.null(hmsc$data$Auto)){
			if(is.null(hmsc$data$X)){
				### Only Random
				for(i in 1:niter){
					PredRandom <- matrix(0, nrow=nsp, ncol=nRandom)
					for(j in 1:nRandom){
						paramLatent<-hmsc$results$estimation$paramLatent[[i,j]]
						nLatentRandom<-ncol(paramLatent)
						
						for(k in 1:nLatentRandom){
							PredRandom[,j] <- PredRandom[,j] + paramLatent[,k]*paramLatent[,k]
						}
					}
				
					variPartRandom <- variPartRandom + PredRandom/replicate(nRandom, rowSums(PredRandom))
				}
				
				variPartRandom <- variPartRandom / niter
				res<-variPartRandom
			}else{
				### X and Random
				for(i in 1:niter){
					predXTotal <- rep(0, nsp)
					predXSplit <- matrix(0, nrow=nsp, ncol=nGroups)
					PredRandom <- matrix(0, nrow=nsp, ncol=nRandom)
					
					if(any(names(hmsc$data) == "Tr")){
						### Calculate R2 related to traits
						predX <- tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
						
						meansParamX <- t(hmsc$results$estimation$paramTr[,,i]%*%hmsc$data$Tr)
						predTr <- tcrossprod(hmsc$data$X,meansParamX)
						traitR2Base <- vector(length=2,"numeric")
						
						for(j in 1:nsite){
							covPred <- cov(cbind(predTr[j,], predX[j,]))
							traitR2Base[1] <- traitR2Base[1] + covPred[1,2]*covPred[2,1]
							traitR2Base[2] <- traitR2Base[2] + covPred[1,1]*covPred[2,2]
						}
						
						traitR2 <- traitR2+traitR2Base[1]/traitR2Base[2]
					}
					
					
					### Calculate R2 related to X
					for(j in 1:nsp){
						predXTotalSub <- hmsc$results$estimation$paramX[j,,i] %*% crossprod(covX,hmsc$results$estimation$paramX[j,,i])
						predXTotal[j] <- predXTotal[j] + predXTotalSub
						for(k in 1:nGroups){
							sel <- groupX==groupXLev[k]
							predXPart <- hmsc$results$estimation$paramX[j,sel,i] %*% crossprod(covX[sel,sel],hmsc$results$estimation$paramX[j,sel,i])
							predXSplit[j,k] <- predXSplit[j,k] + predXPart
						}
					}
					
					### Calculate R2 related to Random
					for(j in 1:nRandom){
						paramLatent<-hmsc$results$estimation$paramLatent[[i,j]]
						nLatentRandom<-ncol(paramLatent)
						
						for(k in 1:nLatentRandom){
							PredRandom[,j] <- PredRandom[,j] + paramLatent[,k]*paramLatent[,k]
						}
					}
				
					variTotal <- predXTotal + rowSums(PredRandom)
					variPartX <- variPartX + predXTotal/variTotal;
					variPartRandom <- variPartRandom + PredRandom/replicate(nRandom, variTotal)
					variPartXSplit <- variPartXSplit + predXSplit/replicate(nGroups, apply(predXSplit, 1, sum));
				}
				
				
				variPartX <- variPartX / niter
				variPartRandom <- variPartRandom / niter
				variPartXSplit <- variPartXSplit / niter
				variPart<-cbind(replicate(nGroups,variPartX)*variPartXSplit,variPartRandom)

				if(any(names(hmsc$data) == "Tr")){
					traitR2 <- traitR2 / niter
					res<-list(variPart=variPart,
							  traitR2=traitR2)
				}else{
					res<-variPart
				}
			}
		}else{
			if(is.null(hmsc$data$X)){
				### Auto and Random
				for(i in 1:niter){
					PredAuto <- matrix(0, nrow=nsp, ncol=nAuto)
					PredRandom <- matrix(0, nrow=nsp, ncol=nRandom)
					### R2 Random
					for(j in 1:nRandom){
						paramLatent<-hmsc$results$estimation$paramLatent[[i,j]]
						nLatentRandom<-ncol(paramLatent)
						
						for(k in 1:nLatentRandom){
							PredRandom[,j] <- PredRandom[,j] + paramLatent[,k]*paramLatent[,k]
						}
					}
					
					### R2 Auto
					for(j in 1:nAuto){
						paramLatentAuto<-hmsc$results$estimation$paramLatentAuto[[i,j]]
						nLatentAuto<-ncol(paramLatentAuto)
						
						for(k in 1:nLatentAuto){
							PredAuto[,j] <- PredAuto[,j] + paramLatentAuto[,k]*paramLatentAuto[,k]
						}
					}
				
					variTotal <- rowSums(cbind(PredRandom,PredAuto))
					variPartRandom <- variPartRandom + PredRandom/replicate(nRandom, variTotal)
					variPartAuto <- variPartAuto + PredAuto/replicate(nAuto, variTotal)
				}
				
				variPartRandom <- variPartRandom / niter
				variPartAuto <- variPartAuto / niter
				
				res<-cbind(variPartRandom,variPartAuto)
				
			}else{
				### X, Auto and Random
				for(i in 1:niter){
					predXTotal <- rep(0, nsp)
					predXSplit <- matrix(0, nrow=nsp, ncol=nGroups)
					PredAuto <- matrix(0, nrow=nsp, ncol=nAuto)
					PredRandom <- matrix(0, nrow=nsp, ncol=nRandom)
					
					if(any(names(hmsc$data) == "Tr")){
						### Calculate R2 related to traits
						predX <- tcrossprod(hmsc$data$X,hmsc$results$estimation$paramX[,,i])
						
						meansParamX <- t(hmsc$results$estimation$paramTr[,,i]%*%hmsc$data$Tr)
						predTr <- tcrossprod(hmsc$data$X,meansParamX)
						traitR2Base <- vector(length=2,"numeric")
						
						for(j in 1:nsite){
							covPred <- cov(cbind(predTr[j,], predX[j,]))
							traitR2Base[1] <- traitR2Base[1] + covPred[1,2]*covPred[2,1]
							traitR2Base[2] <- traitR2Base[2] + covPred[1,1]*covPred[2,2]
						}
						
						traitR2 <- traitR2+traitR2Base[1]/traitR2Base[2]
					}
					
					
					### Calculate R2 related to X
					for(j in 1:nsp){
						predXTotalSub <- hmsc$results$estimation$paramX[j,,i] %*% crossprod(covX,hmsc$results$estimation$paramX[j,,i])
						predXTotal[j] <- predXTotal[j] + predXTotalSub
						for(k in 1:nGroups){
							sel <- groupX==groupXLev[k]
							predXPart <- hmsc$results$estimation$paramX[j,sel,i] %*% crossprod(covX[sel,sel],hmsc$results$estimation$paramX[j,sel,i])
							predXSplit[j,k] <- predXSplit[j,k] + predXPart
						}
					}
					
					### Calculate R2 related to Auto
					for(j in 1:nAuto){
						paramLatentAuto<-hmsc$results$estimation$paramLatentAuto[[i,j]]
						nLatentAuto<-ncol(paramLatentAuto)
						
						for(k in 1:nLatentAuto){
							PredAuto[,j] <- PredAuto[,j] + paramLatentAuto[,k]*paramLatentAuto[,k]
						}
					}
				
					### Calculate R2 related to Random
					for(j in 1:nRandom){
						paramLatent<-hmsc$results$estimation$paramLatent[[i,j]]
						nLatentRandom<-ncol(paramLatent)
						
						for(k in 1:nLatentRandom){
							PredRandom[,j] <- PredRandom[,j] + paramLatent[,k]*paramLatent[,k]
						}
					}
				
					variTotal <- predXTotal + cbind(PredRandom,PredAuto)
					variPartX <- variPartX + predXTotal/variTotal;
					variPartAuto <- variPartAuto + PredAuto/replicate(nAuto, variTotal)
					variPartRandom <- variPartRandom + PredRandom/replicate(nRandom, variTotal)
					variPartXSplit <- variPartXSplit + predXSplit/replicate(nGroups, apply(predXSplit, 1, sum));
				}
				
				
				variPartX <- variPartX / niter
				variPartAuto <- variPartAuto / niter
				variPartRandom <- variPartRandom / niter
				variPartXSplit <- variPartXSplit / niter
				variPart<-cbind(replicate(nGroups,variPartX)*variPartXSplit,variPartRandom,variPartAuto)

				if(any(names(hmsc$data) == "Tr")){
					traitR2 <- traitR2 / niter
					res<-list(variPart=variPart,
							  traitR2=traitR2)
				}else{
					res<-variPart
				}
			}
		}
	}
	
	return(res)
}
