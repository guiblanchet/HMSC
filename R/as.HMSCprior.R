#' @rdname as.HMSCdata
#' @export
as.HMSCprior <-
function(data,family="probit",varDistShape=NULL,varDistScale=NULL,
		 varXDf=NULL,varXScaleMat=NULL,meansParamX=NULL,
		 varMeansParamX=NULL,residVar=NULL,paramTr=NULL,varTr=NULL,
		 paramPhylo=NULL,paramAutoDist=NULL,paramAutoWeight=NULL,
		 shrinkOverall=NULL,shrinkSpeed=NULL,shrinkLocal=NULL){

	### Find the data type in the data object
	dataType <- names(data)

	### Remove "Y"
	if(!is.null(data$Y)){
		dataType <- dataType[-which(dataType=="Y")]
	}

	### Number of datatypes
	nDataType <- length(dataType)

	### Number of variables in X
	if(!is.null(data$X)){
		nparamX <- ncol(data$X)
	}

	### Define parameters for Auto
	if(!is.null(data$Auto)){
		if(!is.null(paramAutoDist) & !is.null(paramAutoWeight)){
			if(nrow(paramAutoDist)!=nrow(paramAutoWeight)){
				stop("'paramAutoDist' should have the same number of rows as 'paramAutoWeight'")
			}
			if(ncol(paramAutoDist)!=ncol(paramAutoWeight)){
				stop("'paramAutoDist' should have the same number of columns as 'paramAutoWeight'")
			}
		}

		if(is.null(paramAutoDist) | is.null(paramAutoWeight)){
			print("'paramAutoDist' and 'paramAutoWeight' were defined, if one was not NULL it was overwritten")
		}
	}

	#================================================
	### If there is one type of explanatory variables
	#================================================
	if(nDataType == 1){
		if(dataType == "X"){
			priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
			priorsResidVar <- flatPriorResidVar(residVar)

			priors <- list(meansParamX=priorsX$meansParamX,
						 varMeansParamX=priorsX$varMeansParamX,
						 varXDf=priorsX$varXDf,
						 varXScaleMat=priorsX$varXScaleMat,
						 residVar=priorsResidVar$residVar)

		}else{
			if(all(dataType == "Random")){
				priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
				priorsResidVar <- flatPriorResidVar(residVar)

				priors <- list(residVar=priorsResidVar$residVar,
							 shrinkOverall=priorsRandom$shrinkOverall,
							 shrinkSpeed=priorsRandom$shrinkSpeed,
							 shrinkLocal=priorsRandom$shrinkLocal)
			}else{
				if(all(dataType == "Auto")){
					priorsResidVar <- flatPriorResidVar(residVar)
					priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

					priors <- list(residVar=priorsResidVar$residVar,
								 paramAutoDist=priorsAuto$paramAutoDist,
								 paramAutoWeight=priorsAuto$paramAutoWeight,
								 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
								 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
								 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
				}
			}
		}
	}

	#==================================================
	### If there are two types of explanatory variables
	#==================================================
	if(nDataType == 2){
		if(all(dataType == c("X","Tr"))){
			nTr <- nrow(data$Tr)
			priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
			priorsResidVar <- flatPriorResidVar(residVar)

			priors <- list(paramTr=priorsXTr$paramTr,
						 varTr=priorsXTr$varTr,
						 varXDf=priorsXTr$varXDf,
						 varXScaleMat=priorsXTr$varXScaleMat,
						 residVar=priorsResidVar$residVar)
		}else{
			if(all(dataType == c("X","Random"))){
				priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
				priorsResidVar <- flatPriorResidVar(residVar)
				priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

				priors <- list(meansParamX=priorsX$meansParamX,
							 varMeansParamX=priorsX$varMeansParamX,
							 varXDf=priorsX$varXDf,
							 varXScaleMat=priorsX$varXScaleMat,
							 residVar=priorsResidVar$residVar,
							 shrinkOverall=priorsRandom$shrinkOverall,
							 shrinkSpeed=priorsRandom$shrinkSpeed,
							 shrinkLocal=priorsRandom$shrinkLocal)
			}else{
				if(all(dataType == c("X","Phylo"))){
					priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
					priorsResidVar <- flatPriorResidVar(residVar)
					priorsPhylo <- flatPriorPhylo(paramPhylo)

					### List of priors
					priors <- list(meansParamX=priorsX$meansParamX,
								 varMeansParamX=priorsX$varMeansParamX,
								 varXDf=priorsX$varXDf,
								 varXScaleMat=priorsX$varXScaleMat,
								 residVar=priorsResidVar$residVar,
								 paramPhylo=priorsPhylo)
				}else{
					if(all(dataType == c("X","Auto"))){
						priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
						priorsResidVar <- flatPriorResidVar(residVar)
						priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

						### List of priors
						priors <- list(meansParamX=priorsX$meansParamX,
									 varMeansParamX=priorsX$varMeansParamX,
									 varXDf=priorsX$varXDf,
									 varXScaleMat=priorsX$varXScaleMat,
									 residVar=priorsResidVar$residVar,
									 paramAutoDist=priorsAuto$paramAutoDist,
									 paramAutoWeight=priorsAuto$paramAutoWeight,
									 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
									 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
									 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
					}else{
						if(all(dataType == c("Auto","Random"))){
							priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
							priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
							priorsResidVar <- flatPriorResidVar(residVar)

						### List of priors
						priors <- list(residVar=priorsResidVar$residVar,
									 shrinkOverall=priorsRandom$shrinkOverall,
									 shrinkSpeed=priorsRandom$shrinkSpeed,
									 shrinkLocal=priorsRandom$shrinkLocal,
									 paramAutoDist=priorsAuto$paramAutoDist,
									 paramAutoWeight=priorsAuto$paramAutoWeight,
									 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
									 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
									 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
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
		if(all(dataType == c("X","Tr","Random"))){
			nTr <- nrow(data$Tr)
			priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
			priorsResidVar <- flatPriorResidVar(residVar)
			priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

			priors <- list(paramTr=priorsXTr$paramTr,
						 varTr=priorsXTr$varTr,
						 varXDf=priorsXTr$varXDf,
						 varXScaleMat=priorsXTr$varXScaleMat,
						 residVar=priorsResidVar$residVar,
						 shrinkOverall=priorsRandom$shrinkOverall,
						 shrinkSpeed=priorsRandom$shrinkSpeed,
						 shrinkLocal=priorsRandom$shrinkLocal)
		}else{
			if(all(dataType == c("X","Tr","Phylo"))){
				nTr <- nrow(data$Tr)
				priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
				priorsResidVar <- flatPriorResidVar(residVar)
				priorsPhylo <- flatPriorPhylo(paramPhylo)

				priors <- list(paramTr=priorsXTr$paramTr,
							 varTr=priorsXTr$varTr,
							 varXDf=priorsXTr$varXDf,
							 varXScaleMat=priorsXTr$varXScaleMat,
							 residVar=priorsResidVar$residVar,
							 paramPhylo=priorsPhylo)
			}else{
				if(all(dataType == c("X","Phylo","Random"))){
					priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
					priorsResidVar <- flatPriorResidVar(residVar)
					priorsPhylo <- flatPriorPhylo(paramPhylo)
					priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

					priors <- list(meansParamX=priorsX$meansParamX,
								 varMeansParamX=priorsX$varMeansParamX,
								 varXDf=priorsX$varXDf,
								 varXScaleMat=priorsX$varXScaleMat,
								 residVar=priorsResidVar$residVar,
								 paramPhylo=priorsPhylo,
								 shrinkOverall=priorsRandom$shrinkOverall,
								 shrinkSpeed=priorsRandom$shrinkSpeed,
								 shrinkLocal=priorsRandom$shrinkLocal)
				}else{
					if(all(dataType == c("X","Tr","Auto"))){
						nTr <- nrow(data$Tr)
						priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
						priorsResidVar <- flatPriorResidVar(residVar)
						priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

						priors <- list(paramTr=priorsXTr$paramTr,
									 varTr=priorsXTr$varTr,
									 varXDf=priorsXTr$varXDf,
									 varXScaleMat=priorsXTr$varXScaleMat,
									 residVar=priorsResidVar$residVar,
									 paramAutoDist=priorsAuto$paramAutoDist,
									 paramAutoWeight=priorsAuto$paramAutoWeight,
									 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
									 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
									 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
					}else{
						if(all(dataType == c("X","Phylo","Auto"))){
							priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
							priorsResidVar <- flatPriorResidVar(residVar)
							priorsPhylo <- flatPriorPhylo(paramPhylo)
							priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

							priors <- list(meansParamX=priorsX$meansParamX,
										 varMeansParamX=priorsX$varMeansParamX,
										 varXDf=priorsX$varXDf,
										 varXScaleMat=priorsX$varXScaleMat,
										 residVar=priorsResidVar$residVar,
										 paramPhylo=priorsPhylo,
										 paramAutoDist=priorsAuto$paramAutoDist,
										 paramAutoWeight=priorsAuto$paramAutoWeight,
										 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
										 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
										 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
						}else{
							if(all(dataType == c("X","Auto","Random"))){
								priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
								priorsResidVar <- flatPriorResidVar(residVar)
								priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
								priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

								### List of priors
								priors <- list(meansParamX=priorsX$meansParamX,
											 varMeansParamX=priorsX$varMeansParamX,
											 varXDf=priorsX$varXDf,
											 varXScaleMat=priorsX$varXScaleMat,
											 residVar=priorsResidVar$residVar,
											 shrinkOverall=priorsRandom$shrinkOverall,
											 shrinkSpeed=priorsRandom$shrinkSpeed,
											 shrinkLocal=priorsRandom$shrinkLocal,
											 paramAutoDist=priorsAuto$paramAutoDist,
											 paramAutoWeight=priorsAuto$paramAutoWeight,
											 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
											 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
											 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
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
		if(all(dataType == c("X","Tr","Phylo","Random"))){
			nTr <- nrow(data$Tr)
			priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
			priorsResidVar <- flatPriorResidVar(residVar)
			priorsPhylo <- flatPriorPhylo(paramPhylo)
			priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

			priors <- list(paramTr=priorsXTr$paramTr,
						 varTr=priorsXTr$varTr,
						 varXDf=priorsXTr$varXDf,
						 varXScaleMat=priorsXTr$varXScaleMat,
						 residVar=priorsResidVar$residVar,
						 paramPhylo=priorsPhylo,
						 shrinkOverall=priorsRandom$shrinkOverall,
						 shrinkSpeed=priorsRandom$shrinkSpeed,
						 shrinkLocal=priorsRandom$shrinkLocal)
		}else{
			if(all(dataType == c("X","Tr","Phylo","Auto"))){
				nTr <- nrow(data$Tr)
				priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
				priorsResidVar <- flatPriorResidVar(residVar)
				priorsPhylo <- flatPriorPhylo(paramPhylo)
				priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

				priors <- list(paramTr=priorsXTr$paramTr,
							 varTr=priorsXTr$varTr,
							 varXDf=priorsXTr$varXDf,
							 varXScaleMat=priorsXTr$varXScaleMat,
							 residVar=priorsResidVar$residVar,
							 paramPhylo=priorsPhylo,
							 paramAutoDist=priorsAuto$paramAutoDist,
							 paramAutoWeight=priorsAuto$paramAutoWeight,
							 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
							 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
							 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
			}else{
				if(all(dataType == c("X","Tr","Auto","Random"))){
					nTr <- nrow(data$Tr)
					priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
					priorsResidVar <- flatPriorResidVar(residVar)
					priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
					priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

					priors <- list(paramTr=priorsXTr$paramTr,
								 varTr=priorsXTr$varTr,
								 varXDf=priorsXTr$varXDf,
								 varXScaleMat=priorsXTr$varXScaleMat,
								 residVar=priorsResidVar$residVar,
								 shrinkOverall=priorsRandom$shrinkOverall,
								 shrinkSpeed=priorsRandom$shrinkSpeed,
								 shrinkLocal=priorsRandom$shrinkLocal,
								 paramAutoDist=priorsAuto$paramAutoDist,
								 paramAutoWeight=priorsAuto$paramAutoWeight,
								 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
								 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
								 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
				}else{
					if(all(dataType == c("X","Phylo","Auto","Random"))){
						priorsX <- flatPriorX(varXDf,varXScaleMat,meansParamX,varMeansParamX,nparamX)
						priorsResidVar <- flatPriorResidVar(residVar)
						priorsPhylo <- flatPriorPhylo(paramPhylo)
						priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
						priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

						### List of priors
						priors <- list(meansParamX=priorsX$meansParamX,
									 varMeansParamX=priorsX$varMeansParamX,
									 varXDf=priorsX$varXDf,
									 varXScaleMat=priorsX$varXScaleMat,
									 residVar=priorsResidVar$residVar,
									 paramPhylo=priorsPhylo,
									 shrinkOverall=priorsRandom$shrinkOverall,
									 shrinkSpeed=priorsRandom$shrinkSpeed,
									 shrinkLocal=priorsRandom$shrinkLocal,
									 paramAutoDist=priorsAuto$paramAutoDist,
									 paramAutoWeight=priorsAuto$paramAutoWeight,
									 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
									 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
									 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
					}
				}
			}
		}
	}

	#===================================================
	### If there are five types of explanatory variables
	#===================================================
	if(nDataType == 5){
		if(all(dataType == c("X","Tr","Phylo","Auto","Random"))){
			nTr <- nrow(data$Tr)
			priorsXTr <- flatPriorXTr(varXDf,varXScaleMat,paramTr,varTr,nTr,nparamX)
			priorsResidVar <- flatPriorResidVar(residVar)
			priorsPhylo <- flatPriorPhylo(paramPhylo)
			priorsRandom <- flatPriorRandom(shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.
			priorsAuto <- flatPriorAuto(data,paramAutoDist,paramAutoWeight,shrinkOverall,shrinkSpeed,shrinkLocal,family) # For now, the argument "family" is not used but later it will be.

			priors <- list(paramTr=priorsXTr$paramTr,
						 varTr=priorsXTr$varTr,
						 varXDf=priorsXTr$varXDf,
						 varXScaleMat=priorsXTr$varXScaleMat,
						 residVar=priorsResidVar$residVar,
						 paramPhylo=priorsPhylo,
						 shrinkOverall=priorsRandom$shrinkOverall,
						 shrinkSpeed=priorsRandom$shrinkSpeed,
						 shrinkLocal=priorsRandom$shrinkLocal,
						 paramAutoDist=priorsAuto$paramAutoDist,
						 paramAutoWeight=priorsAuto$paramAutoWeight,
						 shrinkOverallAuto=priorsAuto$shrinkOverallAuto,
						 shrinkSpeedAuto=priorsAuto$shrinkSpeedAuto,
						 shrinkLocalAuto=priorsAuto$shrinkLocalAuto)
		}
	}

	if(family=="probit"){
		priors <- list(param=priors)
		attributes(priors) <- list(distr="probit")
		names(priors) <- c("param")
	}

	if(family=="poisson"){
		priors <- list(param=priors)
		attributes(priors) <- list(distr="poisson")
		names(priors) <- c("param")
	}

	if(family=="overPoisson"){
		priors <- list(param=priors)
		attributes(priors) <- list(distr="overPoisson")
		names(priors) <- c("param")
	}

	if(family=="gaussian"){
		priors <- list(param=priors)

		attributes(priors) <- list(distr="gaussian")
		names(priors) <- c("param")
	}


	class(priors) <- "HMSCprior"

	return(priors)
}
