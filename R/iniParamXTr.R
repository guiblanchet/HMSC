#' @rdname iniParamX
#' @importFrom stats rnorm
#' @importFrom stats cov
#' @importFrom stats rWishart
#' @importFrom stats coef
#' @importFrom stats glm
#' @importFrom stats binomial
#' @importFrom stats lm
#' @export
iniParamXTr <-
function(data,priors,paramX=NULL,varX=NULL,paramTr=NULL){
	### Number of X variables
	nparamX<-ncol(data$X)
	
	### Number of species
	nsp<-ncol(data$Y)
	#------------------
	### Initiate paramX
	#------------------
	if(nrow(data$Y)!=0){
		if(is.null(paramX)){
			paramX<-matrix(NA,nrow=nsp,ncol=nparamX)
			
			options(warn=-1)
			if(attributes(priors)$distr=="probit"){
				for(i in 1:nsp){
					paramX[i,]<-coef(glm(data$Y[,i]~-1+.,data=as.data.frame(data$X),family=binomial(link = "probit")))
				}
			}
			if(attributes(priors)$distr=="gaussian"){
				for(i in 1:nsp){
					paramX[i,]<-coef(lm(data$Y[,i]~-1+.,data=as.data.frame(data$X)))
				}
			}
			options(warn=0)
			
			### Correct for extreme values
			paramXtoCorrPos<-which(paramX>4,arr.ind=TRUE)
			paramXtoCorrNeg<-which(paramX< -4,arr.ind=TRUE)
			
			if(length(paramXtoCorrPos)>0){
				paramX[paramXtoCorrPos]<- 4
			}
			if(length(paramXtoCorrNeg)>0){
				paramX[paramXtoCorrNeg]<- -4
			}
		}
	}else{
		paramX<-matrix(rnorm(nsp*nparamX),nrow=nsp,ncol=nparamX)
	}
	
	#-------------------
	### Initiate paramTr
	#-------------------
	if(is.null(paramTr)){
		nTr<-nrow(data$Tr)
		paramTr<-matrix(rnorm(nparamX*nTr),nrow=nparamX,ncol=nTr)
	}
	
	#----------------
	### Initiate varX
	#----------------
	if(is.null(varX)){
		varX<-cov(paramX)+diag(0.1,nparamX)
		
		if(any(is.na(varX))){
			precX <- rWishart(1,paramX+1,diag(nparamX))[,,1]
			varX <- solve(precX)
		}else{
			precX<-solve(varX)
		}
	}else{
		precX<-solve(varX)
	}
		
	param<-list(paramX=paramX,
				varX=varX,
				precX=precX,
				paramTr=paramTr)

	### Name objects
	rownames(param$paramX)<-colnames(data$Y)
	colnames(param$paramX)<-colnames(data$X)
	
	rownames(param$varX)<-colnames(data$X)
	colnames(param$varX)<-colnames(data$X)
		
	rownames(param$precX)<-colnames(data$X)
	colnames(param$precX)<-colnames(data$X)
	
	rownames(param$paramTr)<-colnames(data$X)
	colnames(param$paramTr)<-rownames(data$Tr)
	
	return(param)
}
