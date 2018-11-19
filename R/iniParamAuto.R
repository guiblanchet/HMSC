#' @importFrom stats rnorm
#' @importFrom stats rgamma
#' @importFrom stats dist
iniParamAuto <-
function(data,priors,paramAuto=NULL,latentAuto=NULL,paramLatentAuto=NULL,shrinkLocalAuto=NULL,paramShrinkGlobalAuto=NULL){

	### Number of species
	nsp<-ncol(data$Y)

	### Some basic objects about Random
	nAuto<-length(data$Auto) #nr

	### Number of levels in each auto effect considered
	nLevelAuto<-sapply(data$Auto,function(x) nlevels(x[,1]))
	AutoDist<-vector("list",length=nAuto)

	### Calculate distance between groups
	for(i in 1:length(nLevelAuto)){
		nAutoCoord<-ncol(data$Auto[[i]])-1
		AutoCoordMean<-matrix(NA,nrow=nLevelAuto[i],ncol=nAutoCoord)

		for(j in 1:nAutoCoord){
			AutoCoordMean[,j]<-tapply(data$Auto[[i]][,j+1],data$Auto[[i]][,1],mean)
		}

		AutoDist[[i]]<-dist(AutoCoordMean)
	}

	### Find the number of latent variables if some values are given
	if(is.null(latentAuto)){
		nLatentAuto<-max(2,floor(log(nsp)*3))
	}else{
		if(!is.null(paramLatentAuto)){
			nLatentAuto<-sapply(paramLatentAuto,ncol)
		}else{
			if(!is.null(shrinkLocalAuto)){
				nLatentAuto<-sapply(shrinkLocalAuto,ncol)
			}else{
				if(!is.null(paramShrinkGlobalAuto)){
					nLatentAuto<-sapply(paramShrinkGlobalAuto,length)
				}else{
					nLatentAuto<-sapply(latentAuto,ncol)
				}
			}
		}
	}

	#---------------------
	### Initiate paramAuto
	#---------------------
	if(is.null(paramAuto)){
		paramAuto<-vector("list",length=nAuto)
		maxAutoDist<-sapply(AutoDist, max)

		for(i in 1:nAuto){
			paramAuto[[i]][1:nLatentAuto[i]]<-sample(priors$param$paramAutoDist[-1,],nLatentAuto[i])
		}
	}else{
		for(i in 1:nAuto){
			for(j in 1:nLatentAuto[i]){
				paramAuto[[i]][j]<-priors$param$paramAutoDist[which.min(abs(abs(priors$param$paramAutoDist[,1])-paramAuto[[i]][j])),1]
			}
		}
	}

	#----------------------
	### Initiate latentAuto
	#----------------------
	if(is.null(latentAuto)){
		### Construct object to weighted the distance with paramAuto using the exponential function
		wAutoDist<-vector("list",length=nAuto)
		for(i in 1:nAuto){
			wAutoDist[[i]]<-vector("list",length=nLatentAuto[i])
		}

		for(i in 1:nAuto){
			AutoDistMat<-as.matrix(AutoDist[[i]])
			for(j in 1:nLatentAuto[i]){
				wAutoDist[[i]][[j]]<-exp(-AutoDistMat/paramAuto[[i]][j])
			}
		}

		### Object storing latentAuto
		latentAuto<-vector("list",length=nAuto)
		dim(latentAuto)<-c(1,nAuto)
		colnames(latentAuto)<-paste("auto",1:nAuto,sep="")
		for(i in 1:nAuto){
			latentAuto[[1,i]]<-matrix(NA,nrow=nLevelAuto[i],ncol=nLatentAuto)
		}

		for(i in 1:nAuto){
			for(j in 1:nLatentAuto){
				latentAuto[[1,i]][,j]<-rmvnorm(1,rep(0,nrow(wAutoDist[[i]][[j]])),wAutoDist[[i]][[j]])
			}
			rownames(latentAuto[[1,i]])<-levels(data$Auto[[i]][,1])
			colnames(latentAuto[[1,i]])<-paste("latentAuto",1:nLatentAuto[i],sep="")
		}
	}

	#---------------------------
	### Initiate paramLatentAuto
	#---------------------------
	if(is.null(paramLatentAuto)){
		### Object storing paramLatentAuto
		paramLatentAuto<-vector("list",length=nAuto)
		dim(paramLatentAuto)<-c(1,nAuto)
		colnames(paramLatentAuto)<-paste("auto",1:nAuto,sep="")

		for(i in 1:nAuto){
			paramLatentAuto[[1,i]]<-matrix(rnorm(nLatentAuto[i]*nsp),nrow=nsp,ncol=nLatentAuto[i])
			rownames(paramLatentAuto[[1,i]])<-paste("sp",1:nsp,sep="")
			colnames(paramLatentAuto[[1,i]])<-paste("latentAuto",1:nLatentAuto[i],sep="")
		}
	}

	#---------------------------
	### Initiate shrinkLocalAuto
	#---------------------------
	if(is.null(shrinkLocalAuto)){
		shrinkLocalAuto<-vector("list",length=nAuto)
		dim(shrinkLocalAuto)<-c(1,nAuto)
		colnames(shrinkLocalAuto)<-paste("random",1:nAuto,sep="")

		for(i in 1:nAuto){
			shrinkLocalAuto[[1,i]]<-matrix(rgamma(nsp*nLatentAuto[i],priors$param$shrinkLocalAuto/2,priors$param$shrinkLocalAuto/2),nsp,nLatentAuto[i])
			colnames(shrinkLocalAuto[[1,i]])<-paste("latent",1:nLatentAuto[i])
			rownames(shrinkLocalAuto[[1,i]])<-colnames(data$Y)
		}
	}

	#---------------------------------
	### Initiate paramShrinkGlobalAuto
	#---------------------------------
	if(is.null(paramShrinkGlobalAuto)){
		paramShrinkGlobalAuto<-vector("list",length=nAuto)
		dim(paramShrinkGlobalAuto)<-c(1,nAuto)
		colnames(paramShrinkGlobalAuto)<-paste("random",1:nAuto,sep="")

		for(i in 1:nAuto){
			paramShrinkGlobalAuto[[1,i]]<-c(rgamma(1,priors$param$shrinkOverallAuto[1],priors$param$shrinkOverallAuto[2]),rgamma(nLatentAuto[i]-1,priors$param$shrinkSpeedAuto[1],priors$param$shrinkSpeedAuto[2]))
		}
	}

	### List of parameters
	param<-list(paramAuto=paramAuto,
				latentAuto=latentAuto,
				paramLatentAuto=paramLatentAuto,
				shrinkLocalAuto=shrinkLocalAuto,
				paramShrinkGlobalAuto=paramShrinkGlobalAuto)

	### Name object
	names(param$paramAuto)<-names(data$Auto)

	return(param)
}
