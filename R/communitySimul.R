communitySimul<-function(X,Tr=NULL,Random=NULL,nsp=NULL,paramX=NULL,family="probit",paramTr=NULL,meansparamX=NULL,varX=NULL,latent=NULL,paramLatent=NULL,spCor=NULL){
#### F. Guillaume Blanchet - July 2014, October 2014
##########################################################################################
	### This makes the outlierSp inactive
	outlierSp<-NULL
	outlierWeight<-0.001
	#==================================
	### Standard format of the X matrix
	#==================================
	X<-as.matrix(X)
	nsite<-nrow(X)
	nexp<-ncol(X)
	
	if(is.null(colnames(X))){
		colnames(X)<-paste("x",1:nexp,sep="")
	}
	if(is.null(rownames(X))){
		rownames(X)<-paste("site",1:nsite,sep="")
	}
	
	#=====================
	### If paramX is given
	#=====================
	if(!is.null(paramX)){
		nsp<-nrow(paramX)
	}

	#======================
	### If paramTr is given
	#======================
	if(!is.null(paramTr)){
		if(is.null(Tr)){
			stop("'Tr' needs to be given if 'paramTr' is given")
		}
		nsp<-ncol(Tr)
	}
	
	#==========================
	### If paramLatent is given
	#==========================
	if(!is.null(paramLatent)){
		if(is.null(latent) && is.null(Random)){
			stop("If 'paramlatent' is given, 'latent' and 'Random' also needs to be given")
		}
		nsp<-nrow(paramLatent[[1]])
	}
	
	#=====================
	### If latent is given
	#=====================
	if(!is.null(latent)){
		if(is.null(paramLatent) && is.null(Random)){
			stop("If 'latent' is given, 'paramLatent' and 'Random' also needs to be given")
		}
		nsp<-nrow(paramLatent[[1]])
	}
	
	#====================
	### If spCor is given
	#====================
	if(!is.null(spCor)){
		if(!is.null(latent) && !is.null(paramLatent)){
			stop("If 'spCor' is given, 'latent' and 'paramLatent' do not need to be given")
		}
		if(!is.null(Random)){
			stop("If 'spCor' is given, 'Random' does not need to be given")
		}
		nsp<-nrow(spCor)
		
		### Basic check
		if(min(spCor) < -1 | max(spCor) > 1){
			stop("'spCor' is a correlation matrix, so that values in it should range from -1 to 1")
		}
		
		if(nrow(spCor)!=ncol(spCor)){
			stop("'spCor' needs to be a correlation matrix with as many rows as columns")
		}
		
		if(!is.matrix(spCor)){
			spCor<-as.matrix(spCor)
		}
		
		### Add row and column names to spCor
		if(is.null(colnames(spCor))){
			colnames(spCor)<-paste("y",1:nsp,sep="")
		}
		
		if(is.null(rownames(spCor))){
			rownames(spCor)<-paste("y",1:nsp,sep="")
		}
	}
	
	#===================================
	### Standard format of the Tr matrix
	#===================================
	if(is.null(Tr)){
		if(is.null(nsp)){
			stop("Define 'nsp', the number of species to be simulated")
		}
		Tr<-matrix(1,nrow=1,ncol=nsp)
	}else{
		if(!is.vector(Tr)){
			Tr<-as.matrix(Tr)
		}else{
			Tr<-matrix(Tr,nrow=1)
		}
		### Define nsp if it is not yet defined
		if(!is.null(nsp)){
			if(ncol(Tr)!=nsp){
				stop("'nsp' should be the same as the number of column of Tr")
			}
		}else{
			nsp<-ncol(Tr)
		}
	}
	nTr<-nrow(Tr)
	
	if(is.null(colnames(Tr))){
		colnames(Tr)<-paste("y",1:nsp,sep="")
	}
	if(is.null(rownames(Tr))){
		rownames(Tr)<-paste("t",1:nTr,sep="")
	}
	
	#=======================================
	### Standard format of the Random effect
	#=======================================
	if(!is.null(Random)){
		if(is.factor(Random)){
			Random<-data.frame(random1=Random)
			rownames(Random)<-paste("site",1:nsite,sep="")
			### Number or random effects
			nRandom<-1
		}else{
			if(is.data.frame(Random)){
				if(!all(mapply(is.factor,Random))){
					stop("If 'Random' is a data.frame, it should only include factors")
				}
				### Number or random effects
				nRandom<-ncol(Random)

				colnames(Random)<-paste("random",1:nRandom,sep="")
				rownames(Random)<-paste("site",1:nsite,sep="")
			}else{
				stop("'Random' should be a factor or a data.frame")
			}
		}
	}
	
	#========================================
	### Set parameters for Tr and meansparamX
	#========================================
	if(is.null(paramTr)){
		if(is.null(meansparamX)){
			paramTr<-matrix(rnorm(nTr*nexp),nrow=nexp,ncol=nTr)
			rownames(paramTr)<-paste("p",1:nexp,sep="")
			colnames(paramTr)<-paste("t",1:nTr,sep="")
			
			### Community Mean (projector of traits and parameters)
			meansparamX<-matrix(NA,ncol=nexp,nrow=nsp)
			for(i in 1:nsp){
				meansparamX[i,]<-tcrossprod(Tr[,i],paramTr)
			}
			rownames(meansparamX)<-paste("sp",1:nsp,sep="")
			colnames(meansparamX)<-paste("p",1:nexp,sep="")
		}else{
			paramTrAll<-array(dim=c(nexp,nTr,nsp))
			for(i in 1:nsp){
				paramTrAll[,,i]<-meansparamX[i,]%*%matrix(Tr[,i],nrow=1) # matrix() is to account for the case where there is only 1 trait
			}
			### Approximation of paramTr
			paramTr<-apply(paramTrAll,1:2,mean)
		}
	}else{
		if(is.null(meansparamX)){
			### Community Mean (projector of traits and parameters)
			meansparamX<-matrix(NA,ncol=nexp,nrow=nsp)
			for(i in 1:nsp){
				meansparamX[i,]<-tcrossprod(Tr[,i],paramTr)
			}
			rownames(meansparamX)<-paste("sp",1:nsp,sep="")
			colnames(meansparamX)<-paste("p",1:nexp,sep="")
		}else{
			### Check if meansparamX and meansParamTr (means calculated using Tr and paramTr) match
			meansParamTr<-matrix(NA,ncol=nexp,nrow=nsp)
			for(i in 1:nsp){
				meansParamTr[i,]<-tcrossprod(Tr[,i],paramTr)
			}
			
			if(abs(sum(meansparamX-meansParamTr)) > 10^(-4)){
				warnings("There might be a mismatch between 'meansparamX' and the 'means' calculated using 'Tr' and 'paramTr'")
				warning("The values in 'meansparamX' will be used for the remaining data simulation")
			}
			
			rownames(meansparamX)<-paste("sp",1:nsp,sep="")
			colnames(meansparamX)<-paste("p",1:nexp,sep="")
		}
	}
	
	#=========================
	### Define outlier species
	#=========================
	outlierWeightVec<-rep(1,nsp)
	if(!is.null(outlierSp)){
		noutlierSp<-round(nsp*outlierSp)
		outlierWeightVec[1:noutlierSp]<-outlierWeight
	}
	
	#=======================
	### Set parameters for X
	#=======================
	### Mean of the community paramX
	if(is.null(paramX)){
		### Sigma matrix of the community
		if(is.null(varX)){
			varX<-chol2inv(chol(rWishart(1,nexp+2,diag(nexp))[,,1]))
		}else{
			if(!isSymmetric(varX)){
				stop("'varX' is not a symmetric matrix")
			}
		}
		if(is.null(colnames(varX))){
			colnames(varX)<-paste("p",1:nexp,sep="")
		}
		if(is.null(rownames(varX))){
			rownames(varX)<-paste("p",1:nexp,sep="")
		}
		### paramX for each species
		paramX<-matrix(NA,nrow=nsp,ncol=nexp)
		for(i in 1:nsp){
			paramX[i,]<-mvrnorm(1,mu=meansparamX[i,],Sigma=(1/outlierWeightVec[i])*varX)
		}
		colnames(paramX)<-paste("p",1:nexp,sep="")
		rownames(paramX)<-paste("sp",1:nsp,sep="")
	}else{
		if(is.null(meansparamX)){
			meansparamX<-colMeans(paramX)
		}
		if(is.null(varX)){
			varX<-crossprod(sweep(paramX,2,colMeans(paramX),FUN="-"))
		}
		
	}
	
	#===================
	### Model estimation
	#===================
	EstModel<-tcrossprod(X,paramX)
	
	#==================================================
	### Set parameters for the latent part of the model
	#==================================================
	if(!is.null(Random)){
		if(is.null(latent)){
			### Object storing the latent variables and the parameters
			latent<-vector("list",length=nRandom)
			names(latent)<-paste("random",1:nRandom,sep="")
		}
		
		if(is.null(paramLatent)){
			paramLatent<-vector("list",length=nRandom)
			names(paramLatent)<-paste("random",1:nRandom,sep="")
		}
		
		latentCov<-array(dim=c(nsp,nsp,nRandom))
		dimnames(latentCov)[[1]]<-paste("y",1:nsp,sep="")
		dimnames(latentCov)[[2]]<-paste("y",1:nsp,sep="")
		dimnames(latentCov)[[3]]<-paste("random",1:nRandom,sep="")
		
		### Number of levels in each random effect considered
		nLevelRandom<-mapply(nlevels,Random)
		
		### Construct the model
		for(i in 1:nRandom){
			### Define the latent variables and their associated parameters
			if(is.null(latent[[i]])){
				if(is.null(paramLatent[[i]])){
					### Define the number of latent variables in the model (it will vary with the number of random effects)
					nlatent<-2*(1:nRandom)
					
					latent[[i]]<-matrix(rnorm(nLevelRandom[i]*nlatent[i]),nrow=nLevelRandom[i],ncol=nlatent[i])
					rownames(latent[[i]])<-paste("lev",1:nLevelRandom[i],sep="")
					colnames(latent[[i]])<-paste("pl",1:nlatent[i],sep="")
					
					paramLatent[[i]]<-matrix(rnorm(nlatent[i]*nsp),nrow=nsp,ncol=nlatent[i])
					rownames(paramLatent[[i]])<-paste("y",1:nsp,sep="")
					colnames(paramLatent[[i]])<-paste("pl",1:nlatent[i],sep="")
				}
			}
			### Add the random effect to the estimated model
			EstModel<-EstModel+tcrossprod(latent[[i]][Random[,i],],paramLatent[[i]])

			### Calculate the variance-covariance matrix from the latent variables
			latentCov[,,i]<-tcrossprod(paramLatent[[i]])
			
		}
	}
	
	if(!is.null(spCor)){
		EstModel<-EstModel+mvrnorm(n=nsite,mu=rep(0,nsp),Sigma=spCor)
	}
	
	#======================================
	### Construct species occurrence matrix
	#======================================
	
	if(family="probit"){
		Ylatent<-matrix(rnorm(nsite*nsp,mean=as.vector(EstModel),sd=1),nrow=nsite,ncol=nsp)
		Y<-Ylatent
		Y[Y>0]<-1 # * much faster than ifelse()
		Y[Y<0]<-0
	}
	if(family="poisson"){
		Y<-matrix(rpois(nsite*nsp,lambda=as.vector(exp(EstModel))),nrow=nsite,ncol=nsp)
	}
	
	rownames(Y)<-paste("site",1:nsite,sep="")
	colnames(Y)<-paste("y",1:nsp,sep="")
	
	#=============================
	### Return the results objects
	#=============================
	### Data
	data<-list(Y=as.data.frame(Y),X=as.data.frame(X),Tr=as.data.frame(Tr),Random=Random)
	attributes(data)<-list(names=c("Y","X","Tr","Random"),Ypattern="sp")
	
	class(data)<-"HMSCdata"
	
	### Parameters
	if(!is.null(Random)){
		allparam<-list(paramX=paramX,paramTr=paramTr,meansparamX=meansparamX,varX=varX,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov)
	}else{
		if(!is.null(spCor)){
			allparam<-list(paramX=paramX,paramTr=paramTr,meansparamX=meansparamX,varX=varX,outlier=outlierWeightVec,spCor=spCor)
		}else{
			allparam<-list(paramX=paramX,paramTr=paramTr,meansparamX=meansparamX,varX=varX,outlier=outlierWeightVec)
		}
	}
	
	class(allparam)<-"HMSCparam"
	
	### Final result object
	res<-list(data=data,param=allparam,probMat=pnorm(EstModel))
	class(res)<-"communitySimul"
	
	return(res)
}
