#' @importFrom stats rnorm
#' @importFrom stats rgamma
iniParamRandom <-
function(data,priors,latent=NULL,paramLatent=NULL,shrinkLocal=NULL,paramShrinkGlobal=NULL){

	### Number of species
	nsp<-ncol(data$Y)
	nRandom<-ncol(data$Random) #nr

	if(is.null(latent)){
		nLatent<-rep(max(2,floor(log(nsp)*3)),nRandom)
	}else{
		if(!is.null(paramLatent)){
			nLatent<-sapply(paramLatent,ncol)
		}else{
			if(!is.null(shrinkLocal)){
				nLatent<-sapply(shrinkLocal,ncol)
			}else{
				if(!is.null(paramShrinkGlobal)){
					nLatent<-sapply(paramShrinkGlobal,length)
				}
			}
		}
	}

	#----------------------------
	### Initiate latent variables
	#----------------------------
	if(is.null(latent)){
		### Some basic objects about Random
		nRandomLev<-mapply(nlevels,data$Random)
		Random<-sapply(data$Random,as.numeric)

		### Initiate latent
		latent<-vector("list",length=nRandom)
		dim(latent)<-c(1,nRandom)
		colnames(latent)<-paste("random",1:nRandom,sep="")

		for(i in 1:nRandom){
			latent[[1,i]]<-matrix(rnorm(nRandomLev[i]*nLatent[i]),nrow=nRandomLev[i],ncol=nLatent[i])
			colnames(latent[[1,i]])<-paste("latent",1:nLatent[i])
			rownames(latent[[1,i]])<-levels(data$Random[,i])
		}
	}

	#-----------------------
	### Initiate paramLatent
	#-----------------------
	if(is.null(paramLatent)){
		paramLatent<-vector("list",length=nRandom)
		dim(paramLatent)<-c(1,nRandom)
		colnames(paramLatent)<-paste("random",1:nRandom,sep="")

		for(i in 1:nRandom){
			paramLatent[[1,i]]<-matrix(rnorm(nsp*nLatent[i]),nrow=nsp,ncol=nLatent[i])
			colnames(paramLatent[[1,i]])<-paste("latent",1:nLatent[i])
			rownames(paramLatent[[1,i]])<-colnames(data$Y)
		}
	}

	#-----------------------
	### Initiate shrinkLocal
	#-----------------------
	if(is.null(shrinkLocal)){
		shrinkLocal<-vector("list",length=nRandom)
		dim(shrinkLocal)<-c(1,nRandom)
		colnames(shrinkLocal)<-paste("random",1:nRandom,sep="")

		for(i in 1:nRandom){
			shrinkLocal[[1,i]]<-matrix(rgamma(nsp*nLatent[i],priors$param$shrinkLocal/2,priors$param$shrinkLocal/2),nsp,nLatent[i])
			colnames(shrinkLocal[[1,i]])<-paste("latent",1:nLatent[i])
			rownames(shrinkLocal[[1,i]])<-colnames(data$Y)
		}
	}

	#-----------------------------
	### Initiate paramShrinkGlobal
	#------------------------------
	if(is.null(paramShrinkGlobal)){
		paramShrinkGlobal<-vector("list",length=nRandom)
		dim(paramShrinkGlobal)<-c(1,nRandom)
		colnames(paramShrinkGlobal)<-paste("random",1:nRandom,sep="")

		for(i in 1:nRandom){
			paramShrinkGlobal[[1,i]]<-c(rgamma(1,priors$param$shrinkOverall[1],priors$param$shrinkOverall[2]),rgamma(nLatent[i]-1,priors$param$shrinkSpeed[1],priors$param$shrinkSpeed[2]))
		}
	}

	param<-list(latent=latent,
				paramLatent=paramLatent,
				shrinkLocal=shrinkLocal,
				paramShrinkGlobal=paramShrinkGlobal)

	return(param)
}
