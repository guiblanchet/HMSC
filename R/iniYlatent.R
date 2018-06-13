#' @importFrom stats runif
#' @importFrom stats pnorm
#' @importFrom stats qnorm
iniYlatent <-
function(data,param,family,ncount){

	### A few basic objects
	if(!is.null(data$X)){
		nparamX<-ncol(data$X)
	}
	nsp<-ncol(data$Y)
	nsite<-nrow(data$Y)

	#==========================================
	### Initiate a latent Y for gaussian models
	#==========================================
	if(family=="gaussian"){
		if(!any(is.na(data$Y))){
			### Transform Y into a matrix
			Y<-as.matrix(data$Y)

			### Initial latent Y
			Ylatent<-Y
		}
	}

	#=======================================
	### Initiate a latent Y for other models
	#=======================================
	if((family=="gaussian" & any(is.na(data$Y))) | family!="gaussian"){
		### Transform X and Y into a matrix
		Y<-as.matrix(data$Y)
		if(!is.null(data$X)){
			X<-as.matrix(data$X)
		}

		#---------------------
		### Calculate EstModel
		#---------------------
		if(is.null(data$X)){
			if(is.null(data$Random)){
				if(is.null(data$Auto)){
					EstModel<-matrix(0,nsite,nsp)
				}else{
					### Only Auto
					EstModel<-matrix(0,nsite,nsp)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
				}
			}else{
				if(is.null(data$Auto)){
					### Only Random
					EstModel<-matrix(0,nsite,nsp)
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}else{
					### Random and Auto
					EstModel<-matrix(0,nsite,nsp)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}
			}
		}else{
			if(is.null(data$Random)){
				if(is.null(data$Auto)){
					### Only X
					EstModel<-tcrossprod(X,param$param$paramX)
				}else{
					### X and Auto
					EstModel<-tcrossprod(X,param$param$paramX)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
				}
			}else{
				if(is.null(data$Auto)){
					### X and Random
					EstModel<-tcrossprod(X,param$param$paramX)
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}else{
					### X, Random and Auto
					EstModel<-tcrossprod(X,param$param$paramX)
					for(i in 1:length(param$param$latentAuto)){
						AutoModel<-tcrossprod(param$param$latentAuto[[1,i]][data$Auto[[i]][,1],],param$param$paramLatentAuto[[1,i]])
						EstModel<-EstModel+AutoModel
					}
					for(i in 1:length(param$param$latent)){
						RandomModel<-tcrossprod(param$param$latent[[1,i]][data$Random[,i],],param$param$paramLatent[[1,i]])
						EstModel<-EstModel+RandomModel
					}
				}
			}
		}

		if(family == "gaussian"){
			Ylatent <- sampleYlatentNormal(Y,matrix(0,nrow=nsite,ncol=nsp),param$param$residVar,nsp,nsite)
		}

		if(family == "probit"){
			Y0<-Y==0
			Y1<-Y==1
			YNA<-is.na(Y)
			Ylatent<-sampleYlatentProbit(Y0,Y1,YNA,matrix(0,nrow=nsite,ncol=nsp),EstModel,param$param$residVar,nsp,nsite)
		}

		if(family == "poisson" | family == "overPoisson"){
			Ylatent<-sampleYlatentPoisson(Y,matrix(0,nrow=nsite,ncol=nsp),EstModel,param$param$residVar,nsp,nsite)
		}

		if(family == "logit"){
			Ylatent <- sampleYlatentBinomialLogit(Y,matrix(0,nrow=nsite,ncol=nsp),EstModel,ncount,nsp,nsite)
		}
	}

	class(Ylatent)<-"HMSCYlatent"
	return(Ylatent)
}
