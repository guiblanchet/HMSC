as.HMSCparam <-
function(HMSCdata,HMSCprior,paramX=NULL,varX=NULL,paramTr=NULL,varTr=NULL,family=binomial(link = "probit")){
#### F. Guillaume Blanchet - September 2014
##########################################################################################
	### Number of species
	nsp<-ncol(HMSCdata$Y)
	#---------
	### paramX
	#---------
	if(is.null(paramX)){
		options(warn=-1)
		paramX<-matrix(NA,nrow=ncol(HMSCdata$Y),ncol=ncol(HMSCdata$X))
		for(i in 1:nsp){
			paramX[i,]<-coef(glm(HMSCdata$Y[,i]~-1+.,data=as.data.frame(HMSCdata$X),family=family))
		}
		options(warn=0)
		
		### Correct for very extreme values
		paramX[paramX > 3]<- 3 # * much faster than ifelse()
		paramX[paramX < -3]<- -3
	}
	rownames(paramX)<-colnames(HMSCdata$Y)
	colnames(paramX)<-colnames(HMSCdata$X)
	
	#----------
	### paramTr
	#----------
	if(is.null(varTr)){
		varTr<-solve(tcrossprod(HMSCdata$Tr)+solve(HMSCprior$varTr)) # Vn
	}
	rownames(varTr)<-rownames(HMSCdata$Tr)
	colnames(varTr)<-rownames(HMSCdata$Tr)
	
	if(is.null(paramTr)){
		paramTr<-t(varTr%*%(HMSCdata$Tr%*%paramX+tcrossprod(solve(HMSCprior$varTr),HMSCprior$paramTr))) # Un
	}
	
	colnames(paramTr)<-rownames(HMSCdata$Tr)
	rownames(paramTr)<-colnames(HMSCdata$X)
	
	#--------------
	### meansparamX
	#--------------
	meansparamX<-t(paramTr%*%HMSCdata$Tr)
	
	colnames(meansparamX)<-colnames(HMSCdata$X)
	rownames(meansparamX)<-colnames(HMSCdata$Y)
	
	#-------
	### varX
	#-------
	if(is.null(varX)){
		varX<-cov(paramX)+diag(0.1,ncol(HMSCdata$X))
		precX<-solve(varX)
	}

	rownames(varX)<-colnames(HMSCdata$X)
	colnames(varX)<-colnames(HMSCdata$X)
	
	rownames(precX)<-colnames(HMSCdata$X)
	colnames(precX)<-colnames(HMSCdata$X)

	param<-list(paramX=paramX,
				varX=varX,
				precX=precX,
				meansparamX=meansparamX,
				paramTr=paramTr,
				varTr=varTr)

	class(param)<-"HMSCparam"
	return(param)
}
