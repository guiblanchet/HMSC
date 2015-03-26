as.HMSCprior <-
function(HMSCdata,varXDf=NULL,varXScaleMat=NULL,paramTr=NULL,varTr=NULL,shrinkOverall=NULL,shrinkSpeed=NULL,shrinkLocal=NULL){
#### F. Guillaume Blanchet - September 2014
##########################################################################################
	### Prior for varX
	if(is.null(varXDf)){
		varXDf<-ncol(HMSCdata$X)+1 # f0
	}
	
	if(is.null(varXScaleMat)){
		varXScaleMat<-diag(ncol(HMSCdata$X)) # G0
	}
	
	### Prior for paramTr
	if(is.null(paramTr)){
		paramTr<-matrix(0,nrow=ncol(HMSCdata$X),ncol=nrow(HMSCdata$Tr)) # U0
	}
	
	if(is.null(varTr)){
		varTr<-diag(nrow(HMSCdata$Tr)) # V0
	}
	
	if(!is.null(HMSCdata$Random)){
		#-----------------------------------------------------------------------------------
		### Conrow(HMSCdata$Tr)ols the overall shrinkage level (ad1 = a_1 in Bhattacharya and Dunson 2011)
		### c(ad1, bd1)
		### ad1 = Shape parameter of a gamma distribution
		### bd1 = Scale parameter of a gamma distribution
		###
		### * ad1 will likely need to be tweaked
		#-----------------------------------------------------------------------------------
		if(is.null(shrinkOverall)){
			shrinkOverall<-c(10,1) # Shape parameter of a gamma distribution        # ad1, bd1
		}
		
		#-----------------------------------------------------------------------------------------------------------------------
		### Conrow(HMSCdata$Tr)ols how fast the shrinkage increases as the factor number increases (ad2 = a_2 in Bhattacharya and Dunson 2011)
		### c(ad2, bd2)
		### ad2 = Shape parameter of a gamma distribution
		### bd2 = Scale parameter of a gamma distribution
		###
		### * ad2 will likely need to be tweaked
		#-----------------------------------------------------------------------------------------------------------------------
		if(is.null(shrinkSpeed)){
			shrinkSpeed<-c(15,1) # Shape parameter of a gamma distribution        # ad2, bd2
		}
		
		### Hyperparameter for the local shrinkage parameter (nu in Bhattacharya and Dunson 2011)
		if(is.null(shrinkLocal)){
			shrinkLocal<-3 # df
		}
	}
	
	### Prior for the outlier species estimation
	### Add only when we figure out how this works
#	if(is.null(outlierSp)){
#		outlierSp<-4 # nu2
#	}
	
	if(!is.null(HMSCdata$Random)){
		priors<-list(paramTr=paramTr,
					 varTr=varTr,
					 varXDf=varXDf,
					 varXScaleMat=varXScaleMat,
					 shrinkOverall=shrinkOverall,
					 shrinkSpeed=shrinkSpeed,
					 shrinkLocal=shrinkLocal)
#					 outlierSp=outlierSp)
	}
	
	if(is.null(HMSCdata$Random)){
		priors<-list(paramTr=paramTr,
					 varTr=varTr,
					 varXDf=varXDf,
					 varXScaleMat=varXScaleMat)
#					 outlierSp=outlierSp)
	}
	
	class(priors)<-"HMSCprior"
	return(priors)
}
