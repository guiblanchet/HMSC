#' @export
hmscTheoryMultiSp <-
function(data,param=NULL,priors=NULL,focal,niter=2000,nburn=1000,thin=1,verbose=TRUE){
#### F. Guillaume Blanchet - May 2016
##########################################################################################
	### General checks
	if(niter < nburn){
		stop("'niter' should be equal or larger than 'burning'")
	}

	### Handle verbose arguments
	verbose<-iniVerbose(verbose,niter)

	### A few basic objects
	nsp<-ncol(data$Y)
	nsite<-nrow(data$Y)

	### Transform each data into a matrix
	Y<-as.matrix(data$Y)
	if(any(names(data)=="X")){
		nparamX<-ncol(data$X)
		X<-as.matrix(data$X)
	}

	if(any(names(data)=="Tr")){
		stop("This function is meant to be used without 'Tr'")
	}

	if(any(names(data)=="Phylo")){
		stop("This function is meant to be used without 'Phylo'")
	}

	#====================================================
	### Initiate prior values if they have not been given
	#====================================================
	if(is.null(priors)){
		priors<-as.HMSCprior(data,family="probit")
	}

	#=================================================================
	### Initiate starting parameter values if they have not been given
	#=================================================================
	if(is.null(param)){
		param<-as.HMSCparam(data,priors)
	}

	### Degrees of freedom for Wishart distribution to update precX and varX (removed from the for loop)
	if(!is.null(param)){
		varXDf<-priors$param$varXDf+nsp
	}

	#====================================================
	### Initiate basic objects to define latent variables
	#====================================================
	if(any(names(data)=="Random")){
		### Some basic objects about Random
		nRandom<-ncol(data$Random) #nr
		nRandomLev<-mapply(nlevels,data$Random) #np
		Random<-sapply(data$Random,as.numeric)-1

		if(nRandom > 1){
			stop("There should be only 1 random effect")
		}

		if(nRandomLev != nsite){
			stop("The number of levels in the random effect should equal the number of samples")
		}

		### Initial number of latent variables
		nLatent<-sapply(param$param$latent,ncol)

		### Parameters for the adaptation when calculating the number and importance of latent variables
		adapt<-c(1,0.0005) # c(b0,b1)

		### redund[1] (prop) : Proportion of redundant elements within factor loadings
		### redund [2] (epsilon) : Proportion of redundant elements within factor loadings
		redund<-c(1,0.001) # c(prop,epsilon)
	}

	#===================================================================
	### Initiate basic objects to define autocorrelated latent variables
	#===================================================================
	if(any(names(data)=="Auto")){
		### Some basic objects about RandomAuto
		nAuto<-length(data$Auto)
		nLevelAuto<-sapply(data$Auto,function(x) nlevels(x[,1]))

		if(nAuto > 1){
			stop("There should be only 1 autocorrelated random effect")
		}

		if(nLevelAuto != nsite){
			stop("The number of levels in the autocorrelated random effect should equal the number of samples")
		}

		### Construct AutoCoord to be used as Auto in the mcmc functions
		AutoCoord<-vector("list",length=nAuto)

		for(i in 1:nAuto){
			nAutoCoord<-ncol(data$Auto[[i]])-1
			AutoCoordMean<-matrix(NA,nrow=nLevelAuto[i],ncol=nAutoCoord)

			for(j in 1:nAutoCoord){
				AutoCoordMean[,j]<-tapply(data$Auto[[i]][,j+1],data$Auto[[i]][,1],mean)
			}

			AutoCoord[[i]]<-AutoCoordMean
		}

		### Construct RandomAuto
		RandomAuto<-vector("list",length=nAuto)

		for(i in 1:nAuto){
			RandomAuto[[i]]<-data$Auto[[i]][,1]
		}

		### Calculate the number of levels in
		nAutoLev<-mapply(nlevels,RandomAuto) #np

		### Reorganize RandomAuto so that it can be used in the mcmc function
		RandomAuto<-sapply(RandomAuto,as.numeric)-1

		### Initial number of latent variables
		nLatentAuto<-sapply(param$param$latentAuto,ncol)

		### Parameters for the adaptation when calculating the number and importance of latent variables
		adapt<-c(1,0.0005) # c(b0,b1)

		### redund[1] (prop) : Proportion of redundant elements within factor loadings
		### redund [2] (epsilon) : Proportion of redundant elements within factor loadings
		redund<-c(1,0.001) # c(prop,epsilon)
	}

	#=================================
	### Initiate a latent Y for probit
	#=================================
	Ylatent<-iniYlatent(data,param,family="probit")

	#======================
	### Construct the model
	#======================
	### Find the data type in the data object
	dataType<-names(data)

	### Remove "Y"
	dataType<-dataType[-which(dataType=="Y")]

	### Number of datatypes
	nDataType<-length(dataType)

	if(nDataType==3){
		### Construct model with Random and Auto
		if(all(dataType %in% c("X","Auto","Random"))){
			diagMat <- array(NA,dim=c(nsite,nsite,nsp))
		 	for(i in 1:nsp){
				diagMat[,,i]<-diag(Y[,i])
			}

			result<-mcmcProbitTheoryMultiSp(focal-1,
											Y,
										  Ylatent,
										  X,
										  AutoCoord,
										  RandomAuto,
										  Random,
										  param$param$paramX,
										  param$param$meansParamX,
										  param$param$precX,
										  param$param$residVar,
										  param$param$latent,
										  param$param$paramLatent,
										  param$param$shrinkLocal,
										  param$param$paramShrinkGlobal,
										  param$param$paramAuto,
										  param$param$latentAuto,
										  param$param$paramLatentAuto,
										  param$param$shrinkLocalAuto,
										  param$param$paramShrinkGlobalAuto,
										  priors$param$meansParamX,
										  priors$param$varXScaleMat,
										  priors$param$varXDf,
										  priors$param$residVar[1],
										  priors$param$residVar[2],
										  priors$param$paramAutoWeight,
										  priors$param$paramAutoDist,
										  priors$param$shrinkLocal,
										  priors$param$shrinkOverall[1],
										  priors$param$shrinkOverall[2],
										  priors$param$shrinkSpeed[1],
										  priors$param$shrinkSpeed[2],
										  adapt,
										  redund,
										  nAuto,
										  nAutoLev,
										  nLatentAuto,
										  nRandom,
										  nRandomLev,
										  nLatent,
										  nsp,
										  nsite,
										  nparamX,
										  nrow(priors$param$paramAutoWeight),
										  niter,
										  nburn,
										  thin,
										  verbose,
											diagMat)
		}
	}

	### Name all parts of the result
	result<-nameResult(data,priors,result,niter,nburn,thin)

	#=================
	### Output results
	#=================
	res<-list(results=result,data=data)

	class(res)<-c("hmsc","probit")

	return(res)
}
