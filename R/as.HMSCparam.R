#' @rdname as.HMSCdata
#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats rWishart
#' @importFrom stats runif
#' @export
as.HMSCparam <-
function(data, priors, varDist=NULL, residVar=NULL, paramX=NULL, 
		 meansParamX=NULL, varX=NULL, paramTr=NULL, paramPhylo=NULL, 
		 paramAuto=NULL, latentAuto=NULL, paramLatentAuto=NULL, 
		 shrinkLocalAuto=NULL, paramShrinkGlobalAuto=NULL, 
		 latent=NULL, paramLatent=NULL, shrinkLocal=NULL, 
		 paramShrinkGlobal=NULL, nsp=NULL){
	
	### A few basic objects
	if(is.null(nsp)){
		if(!is.null(data$Y)){
			nsp <- ncol(data$Y)
		}else{
			stop("'nsp' or 'data$Y' should be given to estimate parameters")
		}
	}
	
	### Construct bogus Y matrix to name objects
	if(is.null(data$Y)){
		data$Y <- matrix(NA, nrow=0, ncol=nsp)
		colnames(data$Y) <- paste("sp", 1:nsp, sep="")
	}
	
	### Find the data type in the data object
	dataType <- names(data)
	
	### Remove "Y"
	if(!is.null(data$Y)){
		dataType <- dataType[-which(dataType=="Y")]
	}
		
	### Number of datatypes
	nDataType <- length(dataType)
	
	if(attr(priors, "distr")=="gaussian"){
		if(is.null(varDist)){
			varDist <- rep(1, nsp)
		}
	}
	
	#================================================
	### If there is one type of explanatory variables
	#================================================
	if(nDataType == 1){
		if(dataType == "X"){
			basicParamResidVar <- iniParamResidVar(data, residVar)
			basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
			
			### List of parameters
			param <- list(paramX=basicParamX$paramX, 
						varX=basicParamX$varX, 
						precX=basicParamX$precX, 
						meansParamX=basicParamX$meansParamX, 
						residVar=basicParamResidVar$residVar)
		}else{
			if(all(dataType == "Random")){
				basicParamResidVar <- iniParamResidVar(data, residVar)
				basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
				
				### List of parameters
				param <- list(latent=basicParamRandom$latent, 
							paramLatent=basicParamRandom$paramLatent, 
							shrinkLocal=basicParamRandom$shrinkLocal, 
							paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
							residVar=basicParamResidVar$residVar)
			}else{
				if(dataType == "Auto"){
					basicParamResidVar <- iniParamResidVar(data, residVar)
					basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
						
					### List of parameters
					param <- list(paramAuto=basicParamAuto$paramAuto, 
								latentAuto=basicParamAuto$latentAuto, 
								paramLatentAuto=basicParamAuto$paramLatentAuto, 
								shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
								paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
								residVar=basicParamResidVar$residVar)
				}
			}
		}
	}
	
	#==================================================
	### If there are two types of explanatory variables
	#==================================================
	if(nDataType == 2){
		if(all(dataType == c("X", "Tr"))){
			basicParamResidVar <- iniParamResidVar(data, residVar)
			basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
			
			### List of parameters
			param <- list(paramX=basicParamXTr$paramX, 
						varX=basicParamXTr$varX, 
						precX=basicParamXTr$precX, 
						paramTr=basicParamXTr$paramTr, 
						residVar=basicParamResidVar$residVar)
		}else{
			if(all(dataType == c("X", "Phylo"))){
				basicParamResidVar <- iniParamResidVar(data, residVar)
				basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
				basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
				
				### List of priors
				param <- list(paramX=basicParamX$paramX, 
							varX=basicParamX$varX, 
							precX=basicParamX$precX, 
							meansParamX=basicParamX$meansParamX, 
							paramPhylo=basicParamPhylo$paramPhylo, 
							residVar=basicParamResidVar$residVar)
			}else{
				if(all(dataType == c("X", "Random"))){
					basicParamResidVar <- iniParamResidVar(data, residVar)
					basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
					basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
					
					### List of parameters
					param <- list(paramX=basicParamX$paramX, 
								varX=basicParamX$varX, 
								precX=basicParamX$precX, 
								meansParamX=basicParamX$meansParamX, 
								latent=basicParamRandom$latent, 
								paramLatent=basicParamRandom$paramLatent, 
								shrinkLocal=basicParamRandom$shrinkLocal, 
								paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
								residVar=basicParamResidVar$residVar)
				}else{
					if(all(dataType[1:2] == c("X", "Auto"))){
						basicParamResidVar <- iniParamResidVar(data, residVar)
						basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
						basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
						
						### List of parameters
						param <- list(paramX=basicParamX$paramX, 
									varX=basicParamX$varX, 
									precX=basicParamX$precX, 
									meansParamX=basicParamX$meansParamX, 
									paramAuto=basicParamAuto$paramAuto, 
									latentAuto=basicParamAuto$latentAuto, 
									paramLatentAuto=basicParamAuto$paramLatentAuto, 
									shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
									paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
									residVar=basicParamResidVar$residVar)
					}else{
						if(all(dataType[1:2] == c("Auto", "Random"))){
							basicParamResidVar <- iniParamResidVar(data, residVar)
							basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
							basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
					
							### List of parameters
							param <- list(latent=basicParamRandom$latent, 
										paramLatent=basicParamRandom$paramLatent, 
										shrinkLocal=basicParamRandom$shrinkLocal, 
										paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
										paramAuto=basicParamAuto$paramAuto, 
										latentAuto=basicParamAuto$latentAuto, 
										paramLatentAuto=basicParamAuto$paramLatentAuto, 
										shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
										paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
										residVar=basicParamResidVar$residVar)
						}
					}
				}
			}
		}
	}
	#====================================================
	### If there are three types of explanatory variables
	#====================================================
	if(nDataType == 3){
		if(all(dataType == c("X", "Tr", "Phylo"))){
			basicParamResidVar <- iniParamResidVar(data, residVar)
			basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
			basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
			
			### List of parameters
			param <- list(paramX=basicParamXTr$paramX, 
						varX=basicParamXTr$varX, 
						precX=basicParamXTr$precX, 
						paramTr=basicParamXTr$paramTr, 
						paramPhylo=basicParamPhylo$paramPhylo, 
						residVar=basicParamResidVar$residVar)
		}else{
			if(all(dataType == c("X", "Tr", "Auto"))){
				basicParamResidVar <- iniParamResidVar(data, residVar)
				basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
				basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
				
				### List of parameters
				param <- list(paramX=basicParamXTr$paramX, 
							varX=basicParamXTr$varX, 
							precX=basicParamXTr$precX, 
							paramTr=basicParamXTr$paramTr, 
							paramAuto=basicParamAuto$paramAuto, 
							latentAuto=basicParamAuto$latentAuto, 
							paramLatentAuto=basicParamAuto$paramLatentAuto, 
							shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
							paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
							residVar=basicParamResidVar$residVar)
			}else{
				if(all(dataType == c("X", "Phylo", "Auto"))){
					basicParamResidVar <- iniParamResidVar(data, residVar)
					basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
					basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
					basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
					
					### List of parameters
					param <- list(paramX=basicParamX$paramX, 
								varX=basicParamX$varX, 
								precX=basicParamX$precX, 
								meansParamX=basicParamX$meansParamX, 
								paramPhylo=basicParamPhylo$paramPhylo, 
								paramAuto=basicParamAuto$paramAuto, 
								latentAuto=basicParamAuto$latentAuto, 
								paramLatentAuto=basicParamAuto$paramLatentAuto, 
								shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
								paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
								residVar=basicParamResidVar$residVar)
				}else{
					if(all(dataType == c("X", "Tr", "Random"))){
						basicParamResidVar <- iniParamResidVar(data, residVar)
						basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
						basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
				
						### List of parameters
						param <- list(paramX=basicParamXTr$paramX, 
									varX=basicParamXTr$varX, 
									precX=basicParamXTr$precX, 
									paramTr=basicParamXTr$paramTr, 
									latent=basicParamRandom$latent, 
									paramLatent=basicParamRandom$paramLatent, 
									shrinkLocal=basicParamRandom$shrinkLocal, 
									paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
									residVar=basicParamResidVar$residVar)
					}else{
						if(all(dataType == c("X", "Phylo", "Random"))){
							basicParamResidVar <- iniParamResidVar(data, residVar)
							basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
							basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
							basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
					
							### List of parameters
							param <- list(paramX=basicParamX$paramX, 
										varX=basicParamX$varX, 
										precX=basicParamX$precX, 
										meansParamX=basicParamX$meansParamX, 
										paramPhylo=basicParamPhylo$paramPhylo, 
										latent=basicParamRandom$latent, 
										paramLatent=basicParamRandom$paramLatent, 
										shrinkLocal=basicParamRandom$shrinkLocal, 
										paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
										residVar=basicParamResidVar$residVar)
						}else{
							if(all(dataType == c("X", "Auto", "Random"))){
								basicParamResidVar <- iniParamResidVar(data, residVar)
								basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
								basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
								basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
							
								### List of parameters
								param <- list(paramX=basicParamX$paramX, 
											varX=basicParamX$varX, 
											precX=basicParamX$precX, 
											meansParamX=basicParamX$meansParamX, 
											latent=basicParamRandom$latent, 
											paramLatent=basicParamRandom$paramLatent, 
											shrinkLocal=basicParamRandom$shrinkLocal, 
											paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
											paramAuto=basicParamAuto$paramAuto, 
											latentAuto=basicParamAuto$latentAuto, 
											paramLatentAuto=basicParamAuto$paramLatentAuto, 
											shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
											paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
											residVar=basicParamResidVar$residVar)
							}
						}
					}
				}
			}
		}
	}

	#===================================================
	### If there are four types of explanatory variables
	#===================================================
	if(nDataType == 4){
		if(all(dataType == c("X", "Tr", "Phylo", "Auto"))){
			basicParamResidVar <- iniParamResidVar(data, residVar)
			basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
			basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
			basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
			
			### List of parameters
			param <- list(paramX=basicParamXTr$paramX, 
						varX=basicParamXTr$varX, 
						precX=basicParamXTr$precX, 
						paramTr=basicParamXTr$paramTr, 
						paramPhylo=basicParamPhylo$paramPhylo, 
						paramAuto=basicParamAuto$paramAuto, 
						latentAuto=basicParamAuto$latentAuto, 
						paramLatentAuto=basicParamAuto$paramLatentAuto, 
						shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
						paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
						residVar=basicParamResidVar$residVar)
		}else{
			if(all(dataType == c("X", "Tr", "Phylo", "Random"))){
				basicParamResidVar <- iniParamResidVar(data, residVar)
				basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
				basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
				basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
				
				### List of parameters
				param <- list(paramX=basicParamXTr$paramX, 
							varX=basicParamXTr$varX, 
							precX=basicParamXTr$precX, 
							paramTr=basicParamXTr$paramTr, 
							paramPhylo=basicParamPhylo$paramPhylo, 
							latent=basicParamRandom$latent, 
							paramLatent=basicParamRandom$paramLatent, 
							shrinkLocal=basicParamRandom$shrinkLocal, 
							paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
							residVar=basicParamResidVar$residVar)
			}else{
				if(all(dataType == c("X", "Tr", "Auto", "Random"))){
					basicParamResidVar <- iniParamResidVar(data, residVar)
					basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
					basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
					basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
				
					### List of parameters
					param <- list(paramX=basicParamXTr$paramX, 
								varX=basicParamXTr$varX, 
								precX=basicParamXTr$precX, 
								paramTr=basicParamXTr$paramTr, 
								latent=basicParamRandom$latent, 
								paramLatent=basicParamRandom$paramLatent, 
								shrinkLocal=basicParamRandom$shrinkLocal, 
								paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
								paramAuto=basicParamAuto$paramAuto, 
								latentAuto=basicParamAuto$latentAuto, 
								paramLatentAuto=basicParamAuto$paramLatentAuto, 
								shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
								paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
								residVar=basicParamResidVar$residVar)
				}else{
					if(all(dataType == c("X", "Phylo", "Auto", "Random"))){
						basicParamResidVar <- iniParamResidVar(data, residVar)
						basicParamX <- iniParamX(data, priors, paramX, meansParamX, varX)
						basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
						basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
						basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
					
						### List of parameters
						param <- list(paramX=basicParamX$paramX, 
									varX=basicParamX$varX, 
									precX=basicParamX$precX, 
									meansParamX=basicParamX$meansParamX, 
									paramPhylo=basicParamPhylo$paramPhylo, 
									latent=basicParamRandom$latent, 
									paramLatent=basicParamRandom$paramLatent, 
									shrinkLocal=basicParamRandom$shrinkLocal, 
									paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
									paramAuto=basicParamAuto$paramAuto, 
									latentAuto=basicParamAuto$latentAuto, 
									paramLatentAuto=basicParamAuto$paramLatentAuto, 
									shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
									paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
									residVar=basicParamResidVar$residVar)
					}
				}
			}
		}
	}
	
	#===================================================
	### If there are five types of explanatory variables
	#===================================================
	if(nDataType == 5){
		if(all(dataType == c("X", "Tr", "Phylo", "Auto", "Random"))){
			basicParamResidVar <- iniParamResidVar(data, residVar)
			basicParamXTr <- iniParamXTr(data, priors, paramX, varX, paramTr)
			basicParamPhylo <- iniParamPhylo(priors, paramPhylo)
			basicParamRandom <- iniParamRandom(data, priors, latent, paramLatent, shrinkLocal, paramShrinkGlobal)
			basicParamAuto <- iniParamAuto(data, priors, paramAuto, latentAuto, paramLatentAuto, shrinkLocalAuto, paramShrinkGlobalAuto)
			
			### List of parameters
			param <- list(paramX=basicParamXTr$paramX, 
						varX=basicParamXTr$varX, 
						precX=basicParamXTr$precX, 
						paramTr=basicParamXTr$paramTr, 
						paramPhylo=basicParamPhylo$paramPhylo, 
						latent=basicParamRandom$latent, 
						paramLatent=basicParamRandom$paramLatent, 
						shrinkLocal=basicParamRandom$shrinkLocal, 
						paramShrinkGlobal=basicParamRandom$paramShrinkGlobal, 
						paramAuto=basicParamAuto$paramAuto, 
						latentAuto=basicParamAuto$latentAuto, 
						paramLatentAuto=basicParamAuto$paramLatentAuto, 
						shrinkLocalAuto=basicParamAuto$shrinkLocalAuto, 
						paramShrinkGlobalAuto=basicParamAuto$paramShrinkGlobalAuto, 
						residVar=basicParamResidVar$residVar)
		}
	}
	
	if(attr(priors, "distr")=="probit"){
		param <- list(param=param)
		attributes(param) <- list(distr="probit")
		names(param) <- c("param")
	}
	
	if(attr(priors, "distr")=="poisson"){
		param <- list(distr=varDist, 
					param=param)
		
		attributes(param) <- list(distr="gaussian")
		names(param) <- c("distr", "param")
	}
	
	class(param) <- "HMSCparam"
	
	return(param)
}
