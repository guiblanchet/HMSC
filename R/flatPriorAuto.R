flatPriorAuto <-
function(data,paramAutoDist=NULL,paramAutoWeight=NULL,
		 shrinkOverallAuto=NULL,shrinkSpeedAuto=NULL,
		 shrinkLocalAuto=NULL,family="probit"){

	### Prior for autocorrelation part of the model
	if(is.null(paramAutoDist) | is.null(paramAutoWeight)){
		nAuto<-length(data$Auto)

		### Construct prior objects
		paramAutoDist<-matrix(NA,nrow=201,ncol=nAuto)
		colnames(paramAutoDist)<-paste("priorAutoDistance",1:nAuto,sep="")

		paramAutoWeight<-matrix(NA,nrow=201,nAuto)
		colnames(paramAutoWeight)<-paste("priorAutoDistance",1:nAuto,sep="")

		### Construct prior weight (more importance is given to a distance of 0)
		one<-rep(1,nrow(paramAutoDist))
		weight<-0.5*one/(sum(one)-1)
		weight[1]<-0.5

		### Define the most extreme points in the sampled data for each levels
		nLevelAuto<-sapply(data$Auto,function(x) nlevels(x[,1]))
		AutoMin<-vector("list",length=nAuto)
		AutoMax<-vector("list",length=nAuto)

		for(i in 1:length(nLevelAuto)){
			nAutoCoord<-ncol(data$Auto[[i]])-1
			AutoCoordMean<-matrix(NA,nrow=nLevelAuto[i],ncol=nAutoCoord)

			for(j in 1:nAutoCoord){
				AutoCoordMean[,j]<-tapply(data$Auto[[i]][,j+1],data$Auto[[i]][,1],mean)
			}

			AutoMin[[i]]<-apply(AutoCoordMean,2,min)
			AutoMax[[i]]<-apply(AutoCoordMean,2,max)
		}

		for(i in 1:nAuto){
			distance<-seq(0,1,by=0.005)*sqrt(sum((AutoMax[[i]]-AutoMin[[i]])^2))
			paramAutoDist[,i]<-distance
			paramAutoWeight[,i]<-weight
		}
	}else{
		if(paramAutoDist<=0){
			stop("values in 'paramAutoDist' must all be postive")
		}
	}

	### Prior associated to the latent variables
	priorsRandom<-flatPriorRandom(shrinkOverallAuto,shrinkSpeedAuto,shrinkLocalAuto,family) # For now, the argument "family" is not used but later it might be.

	print("The priors for the autocorrelated latent variables should be OK for probit models but not necessarily for other models, be careful")

	priors<-list(paramAutoDist=paramAutoDist,
				 paramAutoWeight=paramAutoWeight,
				 shrinkOverallAuto=priorsRandom$shrinkOverall,
				 shrinkSpeedAuto=priorsRandom$shrinkSpeed,
				 shrinkLocalAuto=priorsRandom$shrinkLocal)

	return(priors)
}
