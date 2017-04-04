flatPriorX <- 
function(varXDf=NULL,varXScaleMat=NULL,meansParamX=NULL,varMeansParamX=NULL,nparamX){
	### Prior for varX
	if(is.null(varXDf)){
		varXDf<-nparamX+1
	}
	if(is.null(varXScaleMat)){
		varXScaleMat<-diag(nparamX)
	}

	### Prior for meansParamX
	if(is.null(meansParamX)){
		meansParamX<-as.matrix(rep(0,nparamX))
	}

	if(is.null(varMeansParamX)){
		varMeansParamX<-diag(nparamX)
	}

	### List of priors
	priors<-list(meansParamX=meansParamX,
				 varMeansParamX=varMeansParamX,
				 varXDf=varXDf,
				 varXScaleMat=varXScaleMat)

	return(priors)
}
