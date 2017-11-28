hmsc.OverPoisson <-
function(data,param=NULL,priors=NULL,niter=2000,nburn=1000,thin=1,verbose=TRUE){
#### F. Guillaume Blanchet - March 2017
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
		Tr<-as.matrix(data$Tr)
	}

	if(any(names(data)=="Phylo")){
		Phylo<-as.matrix(data$Phylo)
		### Construct inverse correlation phylogeny matrix
		iPhylo<-cov2cor(chol2inv(chol(Phylo)))
	}

	#====================================================
	### Initiate prior values if they have not been given
	#====================================================
	if(is.null(priors)){
		priors<-as.HMSCprior(data,family="poisson")
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

		Random <- lapply(data$Random,function(x) factor(x, levels=x))
		Random<-sapply(Random,as.numeric)-1

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
		RandomAuto <- lapply(RandomAuto,function(x) factor(x, levels=x))
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
	### Initiate a latent Y for OverPoisson
	#=================================
	Ylatent<-iniYlatent(data,param,family="overPoisson")

	#======================
	### Construct the model
	#======================
	### Find the data type in the data object
	dataType<-names(data)

	### Remove "Y"
	dataType<-dataType[-which(dataType=="Y")]

	### Number of datatypes
	nDataType<-length(dataType)

	if(nDataType==1){
		### Construct model with X
		if(dataType %in% "X"){
			result<-mcmcOverPoissonX(Y,
								Ylatent,
								X,
								param$param$paramX,
								param$param$meansParamX,
								param$param$precX,
								param$param$residVar,
								priors$param$meansParamX,
								priors$param$varXScaleMat,
								priors$param$varXDf,
								priors$param$residVar[1],
								priors$param$residVar[2],
								nsp,
								nsite,
								nparamX,
								niter,
								nburn,
								thin,
								verbose)
		}

		### Construct model with Random
		if(dataType %in% "Random"){
			result<-mcmcOverPoissonLatent(Y,
									 Ylatent,
									 Random,
									 param$param$residVar,
									 param$param$latent,
									 param$param$paramLatent,
									 param$param$shrinkLocal,
									 param$param$paramShrinkGlobal,
									 priors$param$residVar[1],
									 priors$param$residVar[2],
									 priors$param$shrinkLocal,
									 priors$param$shrinkOverall[1],
									 priors$param$shrinkOverall[2],
									 priors$param$shrinkSpeed[1],
									 priors$param$shrinkSpeed[2],
									 adapt,
									 redund,
									 nRandom,
									 nRandomLev,
									 nLatent,
									 nsp,
									 nsite,
									 niter,
									 nburn,
									 thin,
									 verbose)
		}
		### Construct model with Auto
		if(dataType %in% "Auto"){
			result<-mcmcOverPoissonAuto(Y,
								   Ylatent,
								   AutoCoord,
								   RandomAuto,
								   param$param$residVar,
								   param$param$paramAuto,
								   param$param$latentAuto,
								   param$param$paramLatentAuto,
								   param$param$shrinkLocalAuto,
								   param$param$paramShrinkGlobalAuto,
								   priors$param$residVar[1],
								   priors$param$residVar[2],
								   priors$param$paramAutoWeight,
								   priors$param$paramAutoDist,
								   priors$param$shrinkLocalAuto,
								   priors$param$shrinkOverallAuto[1],
								   priors$param$shrinkOverallAuto[2],
								   priors$param$shrinkSpeedAuto[1],
								   priors$param$shrinkSpeedAuto[2],
								   adapt,
								   redund,
								   nAuto,
								   nAutoLev,
								   nLatentAuto,
								   nsp,
								   nsite,
								   nrow(priors$param$paramAutoWeight),
								   niter,
								   nburn,
								   thin,
								   verbose)
		}
	}

	if(nDataType==2){
		### Construct model with X and Tr
		if(all(dataType %in% c("X","Tr"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTr(Y,
								  Ylatent,
								  X,
								  Tr,
								  param$param$paramX,
								  param$param$paramTr,
								  param$param$precX,
								  param$param$residVar,
								  priors$param$paramTr,
								  priors$param$varTr,
								  priors$param$varXScaleMat,
								  priors$param$varXDf,
								  priors$param$residVar[1],
								  priors$param$residVar[2],
								  nsp,
								  nsite,
								  nparamX,
								  nTr,
								  niter,
								  nburn,
								  thin,
								  verbose)
		}

		### Construct model with X and Random
		if(all(dataType %in% c("X","Random"))){
			result<-mcmcOverPoissonXLatent(Y,
									  Ylatent,
									  X,
									  Random,
									  param$param$paramX,
									  param$param$meansParamX,
									  param$param$precX,
									  param$param$residVar,
									  priors$param$meansParamX,
									  priors$param$varXScaleMat,
									  priors$param$varXDf,
									  priors$param$residVar[1],
									  priors$param$residVar[2],
									  priors$param$shrinkLocal,
									  priors$param$shrinkOverall[1],
									  priors$param$shrinkOverall[2],
									  priors$param$shrinkSpeed[1],
									  priors$param$shrinkSpeed[2],
									  adapt,
									  redund,
									  nRandom,
									  nRandomLev,
									  nLatent,
									  nsp,
									  nsite,
									  nparamX,
									  niter,
									  nburn,
									  thin,
									  verbose)
		}

		### Construct model with X and Random
		if(all(dataType %in% c("X","Auto"))){
			result<-mcmcOverPoissonXAuto(Y,
									Ylatent,
									X,
									AutoCoord,
									RandomAuto,
									param$param$paramX,
									param$param$meansParamX,
									param$param$precX,
									param$param$residVar,
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
									priors$param$shrinkLocalAuto,
									priors$param$shrinkOverallAuto[1],
									priors$param$shrinkOverallAuto[2],
									priors$param$shrinkSpeedAuto[1],
									priors$param$shrinkSpeedAuto[2],
									adapt,
									redund,
									nAuto,
									nAutoLev,
									nLatentAuto,
									nsp,
									nsite,
									nparamX,
									nrow(priors$param$paramAutoWeight),
									niter,
									nburn,
									thin,
									verbose)
		}

		### Construct model with X and Phylo
		if(all(dataType %in% c("X","Phylo"))){

			result<-mcmcOverPoissonXPhylo(Y,
									 Ylatent,
									 X,
									 Phylo,
									 iPhylo,
									 param$param$paramX,
									 param$param$meansParamX,
									 param$param$paramPhylo,
									 param$param$precX,
									 param$param$residVar,
									 priors$param$meansParamX,
									 priors$param$varMeansParamX,
									 priors$param$varXScaleMat,
									 priors$param$varXDf,
									 priors$param$residVar[1],
									 priors$param$residVar[2],
									 matrix(priors$param$paramPhylo[,2],ncol=1),
									 priors$param$paramPhylo[,1],
									 nsp,
									 nsite,
									 nparamX,
									 nrow(priors$param$paramPhylo),
									 niter,
									 nburn,
									 thin,
									 verbose)
		}

		### Construct model with Random and Auto
		if(all(dataType %in% c("Auto","Random"))){

			result<-mcmcOverPoissonAutoLatent(Y,
										 Ylatent,
										 AutoCoord,
										 RandomAuto,
										 Random,
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
										 nrow(priors$param$paramAutoWeight),
										 niter,
										 nburn,
										 thin,
										 verbose)
		}
	}

	if(nDataType==3){
		### Construct model with X, Tr and Random
		if(all(dataType %in% c("X","Tr","Random"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrLatent(Y,
										Ylatent,
										X,
										Tr,
										Random,
										param$param$paramX,
										param$param$paramTr,
										param$param$precX,
										param$param$residVar,
										priors$param$paramTr,
										priors$param$varTr,
										priors$param$varXScaleMat,
										priors$param$varXDf,
										priors$param$residVar[1],
										priors$param$residVar[2],
										priors$param$shrinkLocal,
										priors$param$shrinkOverall[1],
										priors$param$shrinkOverall[2],
										priors$param$shrinkSpeed[1],
										priors$param$shrinkSpeed[2],
										adapt,
										redund,
										nRandom,
										nRandomLev,
										nLatent,
										nsp,
										nsite,
										nparamX,
										nTr,
										niter,
										nburn,
										thin,
										verbose)
		}

		### Construct model with X, Tr and Random
		if(all(dataType %in% c("X","Tr","Auto"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrAuto(Y,
									  Ylatent,
									  X,
									  Tr,
									  AutoCoord,
									  RandomAuto,
									  param$param$paramX,
									  param$param$paramTr,
									  param$param$precX,
									  param$param$residVar,
									  param$param$paramAuto,
									  param$param$latentAuto,
									  param$param$paramLatentAuto,
									  param$param$shrinkLocalAuto,
									  param$param$paramShrinkGlobalAuto,
									  priors$param$paramTr,
									  priors$param$varTr,
									  priors$param$varXScaleMat,
									  priors$param$varXDf,
									  priors$param$residVar[1],
									  priors$param$residVar[2],
									  priors$param$paramAutoWeight,
									  priors$param$paramAutoDist,
									  priors$param$shrinkLocalAuto,
									  priors$param$shrinkOverallAuto[1],
									  priors$param$shrinkOverallAuto[2],
									  priors$param$shrinkSpeedAuto[1],
									  priors$param$shrinkSpeedAuto[2],
									  adapt,
									  redund,
									  nAuto,
									  nAutoLev,
									  nLatentAuto,
									  nsp,
									  nsite,
									  nparamX,
									  nTr,
									  nrow(priors$param$paramAutoWeight),
									  niter,
									  nburn,
									  thin,
									  verbose)
		}

		### Construct model with X, Tr and Phylo
		if(all(dataType %in% c("X","Tr","Phylo"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrPhylo(Y,
									   Ylatent,
									   X,
									   Tr,
									   Phylo,
									   iPhylo,
									   param$param$paramX,
									   param$param$paramTr,
									   param$param$paramPhylo,
									   param$param$precX,
									   param$param$residVar,
									   priors$param$paramTr,
									   priors$param$varTr,
									   priors$param$varXScaleMat,
									   priors$param$varXDf,
									   priors$param$residVar[1],
									   priors$param$residVar[2],
									   matrix(priors$param$paramPhylo[,2],ncol=1),
									   priors$param$paramPhylo[,1],
									   nsp,
									   nsite,
									   nparamX,
									   nTr,
									   nrow(priors$param$paramPhylo),
									   niter,
									   nburn,
									   thin,
									   verbose)
		}

		### Construct model with X, Phylo and Random
		if(all(dataType %in% c("X","Phylo","Random"))){

			result<-mcmcOverPoissonXPhyloLatent(Y,
										   Ylatent,
										   X,
										   Phylo,
										   iPhylo,
										   Random,
										   param$param$paramX,
										   param$param$meansParamX,
										   param$param$precX,
										   param$param$paramPhylo,
										   param$param$residVar,
										   priors$param$meansParamX,
										   priors$param$varMeansParamX,
										   priors$param$varXScaleMat,
										   priors$param$varXDf,
										   priors$param$residVar[1],
										   priors$param$residVar[2],
										   matrix(priors$param$paramPhylo[,2],ncol=1),
										   priors$param$paramPhylo[,1],
										   priors$param$shrinkLocal,
										   priors$param$shrinkOverall[1],
										   priors$param$shrinkOverall[2],
										   priors$param$shrinkSpeed[1],
										   priors$param$shrinkSpeed[2],
										   adapt,
										   redund,
										   nRandom,
										   nRandomLev,
										   nLatent,
										   nsp,
										   nsite,
										   nparamX,
										   nrow(priors$param$paramPhylo),
										   niter,
										   nburn,
										   thin,
										   verbose)
		}

		### Construct model with X, Phylo and Auto
		if(all(dataType %in% c("X","Phylo","Auto"))){

			result<-mcmcOverPoissonXPhyloAuto(Y,
										 Ylatent,
										 X,
										 Phylo,
										 iPhylo,
										 AutoCoord,
										 RandomAuto,
										 param$param$paramX,
										 param$param$meansParamX,
										 param$param$precX,
										 param$param$paramPhylo,
										 param$param$residVar,
										 param$param$paramAuto,
										 param$param$latentAuto,
										 param$param$paramLatentAuto,
										 param$param$shrinkLocalAuto,
										 param$param$paramShrinkGlobalAuto,
										 priors$param$meansParamX,
										 priors$param$varMeansParamX,
										 priors$param$varXScaleMat,
										 priors$param$varXDf,
										 priors$param$residVar[1],
										 priors$param$residVar[2],
										 matrix(priors$param$paramPhylo[,2],ncol=1),
										 priors$param$paramPhylo[,1],
										 priors$param$paramAutoWeight,
										 priors$param$paramAutoDist,
										 priors$param$shrinkLocalAuto,
										 priors$param$shrinkOverallAuto[1],
										 priors$param$shrinkOverallAuto[2],
										 priors$param$shrinkSpeedAuto[1],
										 priors$param$shrinkSpeedAuto[2],
										 adapt,
										 redund,
										 nAuto,
										 nAutoLev,
										 nLatentAuto,
										 nsp,
										 nsite,
										 nparamX,
										 nrow(priors$param$paramPhylo),
										 nrow(priors$param$paramAutoWeight),
										 niter,
										 nburn,
										 thin,
										 verbose)
		}

		### Construct model with Random and Auto
		if(all(dataType %in% c("X","Auto","Random"))){

			result<-mcmcOverPoissonXAutoLatent(Y,
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
										  verbose)
		}
	}

	if(nDataType==4){
		### Construct model with X, Tr, Auto and Random
		if(all(dataType %in% c("X","Tr","Auto","Random"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrAutoLatent(Y,
											Ylatent,
											X,
											Tr,
											AutoCoord,
											RandomAuto,
											Random,
											param$param$paramX,
											param$param$paramTr,
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
											priors$param$paramTr,
											priors$param$varTr,
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
											nTr,
											nrow(priors$param$paramAutoWeight),
											niter,
											nburn,
											thin,
											verbose)
		}

		### Construct model with X, Tr, Auto and Random
		if(all(dataType %in% c("X","Phylo","Auto","Random"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXPhyloAutoLatent(Y,
											   Ylatent,
											   X,
											   Phylo,
											   iPhylo,
											   AutoCoord,
											   RandomAuto,
											   Random,
											   param$param$paramX,
											   param$param$meansParamX,
											   param$param$precX,
											   param$param$paramPhylo,
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
											   priors$param$varMeansParamX,
											   priors$param$varXScaleMat,
											   priors$param$varXDf,
											   priors$param$residVar[1],
											   priors$param$residVar[2],
											   matrix(priors$param$paramPhylo[,2],ncol=1),
											   priors$param$paramPhylo[,1],
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
											   nrow(priors$param$paramPhylo),
											   niter,
											   nburn,
											   thin,
											   verbose)
		}

		### Construct model with X, Tr, Phylo and Random
		if(all(dataType %in% c("X","Tr","Phylo","Random"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrPhyloLatent(Y,
											 Ylatent,
											 X,
											 Tr,
											 Phylo,
											 iPhylo,
											 Random,
											 param$param$paramX,
											 param$param$paramTr,
											 param$param$paramPhylo,
											 param$param$precX,
											 param$param$residVar,
											 priors$param$paramTr,
											 priors$param$varTr,
											 priors$param$varXScaleMat,
											 priors$param$varXDf,
											 priors$param$residVar[1],
											 priors$param$residVar[2],
											 matrix(priors$param$paramPhylo[,2],ncol=1),
											 priors$param$paramPhylo[,1],
											 priors$param$shrinkLocal,
											 priors$param$shrinkOverall[1],
											 priors$param$shrinkOverall[2],
											 priors$param$shrinkSpeed[1],
											 priors$param$shrinkSpeed[2],
											 adapt,
											 redund,
											 nRandom,
											 nRandomLev,
											 nLatent,
											 nsp,
											 nsite,
											 nparamX,
											 nTr,
											 nrow(priors$param$paramPhylo),
											 niter,
											 nburn,
											 thin,
											 verbose)
		}

		### Construct model with X, Tr, Phylo and Auto
		if(all(dataType %in% c("X","Tr","Phylo","Auto"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrPhyloAuto(Y,
										   Ylatent,
										   X,
										   Tr,
										   Phylo,
										   iPhylo,
										   AutoCoord,
										   RandomAuto,
										   param$param$paramX,
										   param$param$paramTr,
										   param$param$precX,
										   param$param$paramPhylo,
										   param$param$residVar,
										   param$param$paramAuto,
										   param$param$latentAuto,
										   param$param$paramLatentAuto,
										   param$param$shrinkLocalAuto,
										   param$param$paramShrinkGlobalAuto,
										   priors$param$paramTr,
										   priors$param$varTr,
										   priors$param$varXScaleMat,
										   priors$param$varXDf,
										   priors$param$residVar[1],
										   priors$param$residVar[2],
										   matrix(priors$param$paramPhylo[,2],ncol=1),
										   priors$param$paramPhylo[,1],
										   priors$param$paramAutoWeight,
										   priors$param$paramAutoDist,
										   priors$param$shrinkLocalAuto,
										   priors$param$shrinkOverallAuto[1],
										   priors$param$shrinkOverallAuto[2],
										   priors$param$shrinkSpeedAuto[1],
										   priors$param$shrinkSpeedAuto[2],
										   adapt,
										   redund,
										   nAuto,
										   nAutoLev,
										   nLatentAuto,
										   nsp,
										   nsite,
										   nparamX,
										   nTr,
										   nrow(priors$param$paramPhylo),
										   nrow(priors$param$paramAutoWeight),
										   niter,
										   nburn,
										   thin,
										   verbose)
		}
	}

	if(nDataType==5){
		### Construct model with X, Tr, Phylo, Auto and Random
		if(all(dataType %in% c("X","Tr","Phylo","Auto","Random"))){
			### Basic objects
			nTr<-nrow(data$Tr)

			result<-mcmcOverPoissonXTrPhyloAutoLatent(Y,
												 Ylatent,
												 X,
												 Tr,
												 Phylo,
												 iPhylo,
												 AutoCoord,
												 RandomAuto,
												 Random,
												 param$param$paramX,
												 param$param$paramTr,
												 param$param$precX,
												 param$param$paramPhylo,
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
												 priors$param$paramTr,
												 priors$param$varTr,
												 priors$param$varXScaleMat,
												 priors$param$varXDf,
												 priors$param$residVar[1],
												 priors$param$residVar[2],
												 matrix(priors$param$paramPhylo[,2],ncol=1),
												 priors$param$paramPhylo[,1],
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
												 nTr,
												 nrow(priors$param$paramPhylo),
												 nrow(priors$param$paramAutoWeight),
												 niter,
												 nburn,
												 thin,
												 verbose)
		}
	}

	### Name all parts of the result
	result<-nameResult(data,priors,result,niter,nburn,thin)

	#=================
	### Output results
	#=================
	res<-list(results=result,data=data)

	class(res)<-c("hmsc","overPoisson")

	return(res)
}
