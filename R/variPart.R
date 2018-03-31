#' @title Variation partitioning
#'
#' @description Partition the variation of models
#'
#' @param hmsc An object of the class \code{hmsc}.
#' @param groupX A vector defining how the covariates (\code{X}) are grouped. This argument is ignored if the models does not have any covariates.
#' @param family A character string defining the type of modelling approach to use (See details).
#' @param HMSCprior An object of class HMSCprior. This object is used only for type III variation partitioning analysis, it is not considered otherwise.
#' @param type A character string defining the sum of squares used to perform the variation partitioning. Use "I" for type I sums of squares (default) or "III" for type I sums of squares.
#' @param verbose Logical. Whether comments about the number of submodels to estimate are printed in the console or not. This is only for type III variation partitioning. Default is \code{TRUE}.
#' @param \dots Arguments passed from \code{\link{hmsc}} for the calculation of type III variation partitioning analysis.
#' 
#' @details
#'
#' When deciding how to group the covariates using \code{groupX}, it is essential to also account for the intercept (if present). If it is not included, an error message will be sent. If it is not clear in which group the intercept need to be considered, it is a good idea to include it in a seperate group.
#' 
#' The choice of family is currently limited to "probit", "logit", "gaussian", "poisson" and "overPoisson". This needs to be specified for type III variation partitioning.
#' 
#' If \code{HMSCprior} is \code{NULL} the default flat prior will be used to construct the different sub-models.
#' 
#' Currently, two ways to calculate the variation partitionning have been implemented in this function, type I and type III, which refer to the way the sums of squares are accounted for in the calculation of the variation partitionning. Type I is the approach proposed in Ovaskainen et al. (2017) while type III is the one presented in Leibold et al. (submitted). 
#' 
#' In the calculation of the type III sum of squares, the approach used is similar to the one proposed by Borcard et al. (1992). The full model is the one given in the \code{hmsc} object. To calculate type III variation partitioning, a new set of models will be constructed where the different combinations of datasets will be considered to partition the variation. For this reason, the number of iterations (\code{niter} in \code{\link{hmsc}}), the number of burn-in steps (\code{nburn} in \code{\link{hmsc}}) and the thining (\code{thin} in \code{\link{hmsc}}) for each submodel part will be based on how the full model was constructed. Lastly, the groups of explanatory variables definned by \code{groupX} will be considered here as an independent submodel. So for these reasons, it should be expected that the time it takes to contruct a type III variation partitioning will be much longer than for a type I variation partitioning.
#' 
#' The calculation of type III variation partitioning relies on coefficient of determinations (R^2) to characterize the explained variation in the data. As such, for the moment, it has only been implemented for Gaussian, probit and logistic models.
#' 
#' When calculating type III variation partitioning, the number of iterations, burnin iterations and the thinning is obtained from the \code{hmsc} object. As such, in some circumstance, particularly when there was thining in the original model, it is possible that the number iterations (for estimation and burnin) diverge slightly from the ones used to estimate the parameter of the full model.
#' 
#' @return
#'
#' A list presenting how the explained variation in the data is partitioned between the different groups of the covariates (first part of thet list) and the random effects (second part of the list). Note that both autocorrelated (\code{Auto}) and non-autocorrelated (\code{Random}) are considered here. In addition, this function also calculate the amount of variation explained by traits.
#'
#' @author Guillaume Blanchet (type III), Gleb Tikhonov (type I)
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
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE,  interceptTr = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 2000, nburn = 1000, thin = 1, verbose = 100)
#'
#' #==================================
#' ### Calculate variance partitioning
#' #==================================
#' vp <- variPart(modelDesc, groupX = c("A", "B"), type = "I")
#' vp
#'
#' #===================
#' ### Plot the results
#' #===================
#' barplot(t(vp), las = 2)
#'
#' @keywords univar, multivariate, regression
#' @export
variPart<-function(hmsc, groupX, HMSCprior = NULL, type = "I",...){
  
  ### Basic check
  if(type == "I" | type == "III"){
    stop("type needs to be either 'I' or 'III'")
  }
  
  ### Basic objects
	nsite<-nrow(hmsc$data$Y)
  nsp<-ncol(hmsc$data$Y)
  nX<-ncol(hmsc$data$X)
  nAuto<-length(hmsc$data$Auto)
  nRandom<-ncol(hmsc$data$Random)
  
  ### Check groupX
  if(!is.null(nX)){
  	if(length(groupX) != nX){
  		stop("groupX should have the same length as there are variables in X (including the intercept)")
  	}
  }
  
  #========================
  ### Type I sum of squares
  #========================
  if(type == "I"){
    ### Objects related to X
  	if(!is.null(nX)){
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
  		niter<-nrow(hmsc$results$estimation$paramLatentAuto)
  	}
  	if(!is.null(nRandom)){
  		niter<-nrow(hmsc$results$estimation$paramLatent)
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
  
  					variTotal <- predXTotal + rowSums(cbind(PredRandom,PredAuto))
  					variPartX <- variPartX + predXTotal/variTotal;
  					variPartAuto <- variPartAuto + PredAuto/replicate(nAuto, variTotal)
  					variPartRandom <- variPartRandom + PredRandom/replicate(nRandom, variTotal)
  					variPartXSplit <- variPartXSplit + predXSplit/replicate(nGroups, apply(predXSplit, 1, sum))
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
  }
  
  #==========================
  ### Type III sum of squares
  #==========================
  if(type == "III"){
    ### General checks
    if(is.null(family)){
      stop("'family' needs to be specified")
    }
    
    if(!any(family == c("probit", "logit","gaussian"))){
      stop("family' should be either 'probit', 'logit','gaussian'")
    }
    
    if(family == "logit"){
      if(!all(unique(model$data$Y) == c(0,1))){
        stop("Values should be either 0 or 1 for the logit model")
      }
    }
      
    ### Rename hmsc to prevent confusion
    model <- hmsc
    remove(hmsc)
    
    ### Names "explanatory" variables
    groupXNames <- unique(groupX)
    RandomNames <- names(model$data$Random)
    AutoNames <- names(model$data$Auto)
    
    ### Count submodels to build
    setsVar <- c(groupXNames, RandomNames, AutoNames)
    nsetsVar <- length(setsVar)

    ### Build submodel basis
    subModelBase <- vector("list", length = nsetsVar-1)
    names(subModelBase) <- paste("var",1:(nsetsVar-1),sep="")
    
    for(i in 1:(nsetsVar-1)){
      subModelBase[[i]] <- combn(setsVar,i)
    }
    
    #================================================
    ### Build sets of X matrices to use for modelling
    #================================================
    X <- vector("list", length = nsetsVar-1)
    names(X) <- paste("var",1:(nsetsVar-1),sep="")
    
    for(i in 1:(nsetsVar-1)){
      X[[i]] <- vector("list" , length = ncol(subModelBase[[i]]))
    }
    
    for(i in 1:(nsetsVar-1)){
      for(j in 1:ncol(subModelBase[[i]])){
        X[[i]][[j]] <- matrix(NA,nrow=nrow(model$data$X),ncol=0)
      }
    }
    
    #===================================================
    ### Build sets of Random object to use for modelling
    #===================================================
    Random <- vector("list", length = nsetsVar-1)
    names(Random) <- paste("var",1:(nsetsVar-1),sep="")
    
    for(i in 1:(nsetsVar-1)){
      Random[[i]] <- vector("list" , length = ncol(subModelBase[[i]]))
    }
    
    for(i in 1:(nsetsVar-1)){
      for(j in 1:ncol(subModelBase[[i]])){
        Random[[i]][[j]] <- as.data.frame(matrix(NA,nrow=nrow(model$data$Random),ncol=nRandom))
        colnames(Random[[i]][[j]]) <- RandomNames
      }
    }
    
    #=================================================
    ### Build sets of Auto object to use for modelling
    #=================================================
    Auto <- vector("list", length = nsetsVar-1)
    names(Auto) <- paste("var",1:(nsetsVar-1),sep="")
    
    for(i in 1:(nsetsVar-1)){
      Auto[[i]] <- vector("list" , length = ncol(subModelBase[[i]]))
    }
    
    for(i in 1:(nsetsVar-1)){
      for(j in 1:ncol(subModelBase[[i]])){
        Auto[[i]][[j]] <- vector("list",length=nAuto)
      }
    }
    
    #===============================
    ### Build basis for HMSC objects
    #===============================
    for(i in 1:(nsetsVar-1)){
      for(j in 1:ncol(subModelBase[[i]])){
        ### X
        if(any(subModelBase[[i]][,j] %in% groupXNames)){
          nX <- sum(subModelBase[[i]][,j] %in% groupXNames)
          XSel <- which(groupXNames %in% subModelBase[[i]][,j])
          namesX <- character()
          for(k in 1:nX){
            X[[i]][[j]] <- cbind(X[[i]][[j]],model$data$X[,which(groupX == groupXNames[XSel[k]])])
            namesX <- c(namesX,colnames(model$data$X)[which(groupX == groupXNames[XSel[k]])])
          }
          colnames(X[[i]][[j]]) <- namesX
        }else{
          X[[i]][[j]] <- list(NULL)
        }

        ### Random
        if(any(subModelBase[[i]][,j] %in% RandomNames)){
          nRandom <- sum(subModelBase[[i]][,j] %in% RandomNames)
          RandomSel <- which(RandomNames %in% subModelBase[[i]][,j])
          for(k in 1:nRandom){
            Random[[i]][[j]][,k] <- model$data$Random[,which(RandomNames == RandomNames[RandomSel[k]])]
          }
        }else{
          Random[[i]][[j]] <- list(NULL)
        }
        
        if(any(is.na(Random[[i]][[j]]))){
          colToRm <- unique(which(is.na(Random[[i]][[j]]),arr.ind = TRUE)[,2])
          Random[[i]][[j]] <- Random[[i]][[j]][,-colToRm]
        }
        
        ### Auto
        if(any(subModelBase[[i]][,j] %in% AutoNames)){
          nAuto <- sum(subModelBase[[i]][,j] %in% AutoNames)
          AutoSel <- which(AutoNames %in% subModelBase[[i]][,j])
          for(k in 1:nAuto){
            Auto[[i]][[j]][[k]] <- model$data$Auto[[which(AutoNames == AutoNames[AutoSel[k]])]]
            names(Auto[[i]][[j]])[k] <- names(model$data$Auto)[k]
          }
        }
        
        if(length(Auto[[i]][[j]]) > 1){
          AutoNull <- which(sapply(Auto[[i]][[j]], is.null))
          if(length(AutoNull) > 0){
            Auto[[i]][[j]][AutoNull] <- NULL
          }
        }
      }
    }
    

    #==================================
    ### Calculate niter, nburn and thin
    #==================================
    ### Function to check if a number is a whole number
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
      abs(x - round(x)) < tol
    }
    
    ### From paramX
    if(any(names(model$data) == "X")){
      ### Find total number of iterations
      iterNames <- dimnames(model$results$estimation$paramX)[[3]]
      iterNum <- strsplit(iterNames, "iter")
      iterVal <- as.numeric(sapply(iterNum, function(x) x[2]))
      niterModel <- max(iterVal)
      
      ### Find number of burnin iterations
      burnNames <- dimnames(model$results$burning$paramX)[[3]]
      burnNum <- strsplit(burnNames, "iter")
      burnVal <- as.numeric(sapply(burnNum, function(x) x[2]))
      nburnModel <- max(burnVal)

      ### Find thining
      thinModel <- iterVal[2] - iterVal[1]
    }
    
    ### From paramX
    if(any(names(model$data) == "Random")){
      ### Find total number of iterations
      iterNames <- rownames(model$results$estimation$latent)
      iterNum <- strsplit(iterNames, "iter")
      iterVal <- as.numeric(sapply(iterNum, function(x) x[2]))
      niterModel <- max(iterVal)
      
      ### Find number of burnin iterations
      burnNames <- rownames(model$results$burning$latent)
      burnNum <- strsplit(burnNames, "iter")
      burnVal <- as.numeric(sapply(burnNum, function(x) x[2]))
      nburnModel <- max(burnVal)
      
      ### Find thining
      thinModel <- iterVal[2]-iterVal[1]
    }
    
    if(any(names(model$data) == "Auto")){
      ### Find total number of iterations
      iterNames <- rownames(model$results$estimation$latentAuto)
      iterNum <- strsplit(iterNames, "iter")
      iterVal <- as.numeric(sapply(iterNum, function(x) x[2]))
      niterModel <- max(iterVal)
      
      ### Find number of burnin iterations
      burnNames <- rownames(model$results$burning$latentAuto)
      burnNum <- strsplit(burnNames, "iter")
      burnVal <- as.numeric(sapply(burnNum, function(x) x[2]))
      nburnModel <- max(burnVal)
      
      ### Find thining
      thinModel <- iterVal[2]-iterVal[1]
    }

    ### Check to make sure that niterModel/thinModel is a whole number and correct if necesary
    if(!is.wholenumber(niterModel/thinModel)){
      niterModel <- niterModel + thinModel - 1
    }
    
    if(!is.wholenumber(nburnModel/thinModel)){
      nburnModel <- nburnModel + thinModel - 1
    }
    
    #================
    ### Result object
    #================
    R2model <- vector("list", length = nsetsVar)
    names(R2model) <- paste("var",1:nsetsVar,sep="")
    
    for(i in 1:(nsetsVar-1)){
      R2model[[i]] <- matrix(NA, nrow = nsp, ncol = ncol(subModelBase[[i]]))
      rownames(R2model[[i]]) <- colnames(model$data$Y)
      columnNames <- character()
      for(j in 1:ncol(subModelBase[[i]])){
        columnNames[j] <- paste(subModelBase[[i]][,j],collapse="-")
      }
      colnames(R2model[[i]]) <- columnNames
    }
    
    R2model[[nsetsVar]] <- matrix(NA, nrow = nsp, ncol = 1)
    colnames(R2model[[nsetsVar]]) <- paste(setsVar, collapse = "-")
    rownames(R2model[[nsetsVar]]) <- colnames(model$data$Y)
    
    #=================================================================
    ### Print message stating how many submodels need to be calculated
    #=================================================================
    if(verbose){
      print(paste(sum(sapply(subModelBase,ncol)),"submodels need to be estimated"))
    }
    
    #==========================================================================
    ### Estimate the submodels and calculate the R2 associated to each submodel
    #==========================================================================
    counter <- 1
    for(i in 1:(nsetsVar-1)){
      for(j in 1:ncol(subModelBase[[i]])){
        ### Organize X for HMSCdata
        if(all(sapply(X[[i]][[j]],is.null))){
          XUse <- NULL
        }else{
          XUse <- X[[i]][[j]]
        }

        ### Organize Random for HMSCdata
        if(all(sapply(Random[[i]][[j]],is.null))){
          RandomUse <- NULL
        }else{
          RandomUse <- Random[[i]][[j]]
        }

        ### Organize Auto for HMSCdata
        if(all(sapply(Auto[[i]][[j]],is.null))){
          AutoUse <- NULL
        }else{
          AutoUse <- Auto[[i]][[j]]
        }
        
        ### Make sure no warnings are given
        options(warn = -1)
        
        ### Build HMSCdata object
        HMSCdataObj <- as.HMSCdata(Y = model$data$Y, 
                                   X = XUse, 
                                   Random = RandomUse,
                                   Auto = AutoUse,
                                   interceptX = TRUE,
                                   scaleX = FALSE)
        
        ### Reset warning screen print to default
        options(warn = 0)
        
        ### Estimate parameters for the different submodels
        submodel <- hmsc(data = HMSCdataObj, priors = HMSCprior,
                         family = family, niter = niterModel, 
                         nburn = nburnModel, thin = thinModel, 
                         verbose = FALSE, ...)

        R2model[[i]][,j] <- Rsquared(submodel, adjust = TRUE,  averageSp = FALSE)
        
        if(verbose){
          print(paste("Number of submodels estimated:", counter))
        }
        counter <- counter + 1
      }
    }
    
    ### Calculate R2 for full model
    R2model[[nsetsVar]][,1] <- Rsquared(model,averageSp = FALSE)
    
    #======================
    ### Calculate fractions
    #======================
    ### Build result object
    fraction <- vector("list", length = nsetsVar-1)
    names(fraction) <- paste("overlap",(nsetsVar-1):1,sep="")
    for(i in 1:(nsetsVar-1)){
      fraction[[i]] <- matrix(NA, nrow = nsp, ncol = ncol(subModelBase[[i]]))
      rownames(fraction[[i]]) <- colnames(model$data$Y)
    }
    
    ### Calculate independent fraction (step 1)
    ref <- unlist(strsplit(colnames(R2model[[nsetsVar]]),"-"))
    for(i in 1:(nsetsVar-1)){
      columnNames <- character()
      for(j in 1:ncol(subModelBase[[i]])){
        comp <- unlist(strsplit(colnames(R2model[[i]])[j],"-"))
        fraction[[i]][,j] <- R2model[[nsetsVar]] -  R2model[[i]][,j]
        columnNames[j] <- paste(ref[which(!(ref %in% comp))], collapse = "-")
      }
      colnames(fraction[[i]]) <- columnNames
    }
    
    ### Isolate independent fractions
    for(i in 1:(nsetsVar-2)){
      for(j in 1:ncol(fraction[[i]])){
        fracToSub <- unlist(strsplit(colnames(fraction[[i]])[j], "-"))
        fracSel <- colnames(fraction[[nsetsVar-1]]) %in% fracToSub
        fracSub <- fraction[[nsetsVar-1]][,fracSel]
        if(is.null(ncol(fracSub))){
          fracSub <- matrix(fracSub, nrow = ncol(model$data$Y))
        }
        fraction[[i]][,j] <- fraction[[i]][,j] - rowSums(fracSub)
      }
    }
    
    ### Calculate final fractio (the one that overlaps them all!)
    overlapAll <- drop(R2model[[nsetsVar]]) - rowSums(sapply(fraction,rowSums))
    overlapAll <- matrix(overlapAll,ncol=1)
    rownames(overlapAll) <- colnames(model$data$Y)
    colnames(overlapAll) <- paste(setsVar, collapse = "-")
    
    fraction[[paste("overlap",nsetsVar,sep="")]] <- overlapAll
    fraction <- fraction[order(names(fraction))]
    
    ### Final result
    res <- fraction
  }

  return(res)
}

