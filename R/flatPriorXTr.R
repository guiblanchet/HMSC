flatPriorXTr <-
function(varXDf=NULL,varXScaleMat=NULL,
		 paramTr=NULL,varTr=NULL,nTr,nparamX){

	### Prior for varX
	if(is.null(varXDf)){
		varXDf<-nparamX+1
	}
	if(is.null(varXScaleMat)){
		varXScaleMat<-diag(nparamX)
	}

	### Prior for paramTr
	if(is.null(paramTr)){
		paramTr<-matrix(0,nrow=1,ncol=nTr*nparamX)
	}
	if(is.null(varTr)){
		varparamTr<-diag(nTr*nparamX)
	}

	priors<-list(paramTr=paramTr,
				 varTr=varparamTr,
				 varXDf=varXDf,
				 varXScaleMat=varXScaleMat)

	return(priors)
}
