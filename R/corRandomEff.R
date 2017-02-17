#' @title Species by species  covariance (or correlation) matrix from the parameters of the latent variables (\code{paramLatent})
#'
#' @description Construct a species by species covariance (or correlation) matrix from the parameters of the latent variables (\code{paramLatent}). If multiple sets of latent variables were estimated, this function will calculate as many species by species covariance (or correlation) matrices as the number of sets of latent variables estimated in the \code{hmsc} object. 
#'
#' @param hmsc Object of class hmsc.
#' @param cor Logical. Whether correlation (\code{TRUE}) or covariance (\code{FALSE}) matrices should be constructed. Default is \code{TRUE}.
#' @param latentToUse A vector of integer defining which set of parameters of the estimated non-autocorrelated latent variables should be considered when constructing the correlation matrix. Default uses all the estimated parameters.
#' @param latentAutoToUse A vector of integer defining which set of parameters of the estimated autocorrelated latent variables should be considered when constructing the correlation matrix. Default uses all the estimated parameters.
#' @param burning Logical. Whether the burning iterations should be include (\code{TRUE}) or not (\code{FALSE}).
#'
#' @details
#'
#' The output of this function is an array.
#'
#' @return
#'
#' An object of class \code{corRandomEff}. This object is an array of four dimensions. The correlation between two species (first and second dimension) can be found for a particular iteration (third dimension) associated to a specific random effect.
#'
#' @author F. Guillaume Blanchet
#'
#' @seealso \code{\link{hmsc}}
#'
#' @examples
#'
#' #==================
#' ### Simulating data
#' #==================
#' desc <- matrix(1, nrow = 50, ncol = 1)
#' random1 <- as.factor(1:50)
#' random2 <- as.factor(rep(letters[1:2], each = 25))
#' randEff <- data.frame(rand1 = random1, rand2 = random2)
#' nspecies <- 10
#' 
#' commRandEff2 <- communitySimul(X = desc, Random = randEff, nsp = nspecies)
#' 
#' #=================================
#' ### Formatting data and parameters
#' #=================================
#' formdata <- as.HMSCdata(Y = commRandEff2$data$Y, X = desc, Random = randEff, 
#' 						   interceptX = FALSE)
#' 
#' #==============
#' ### Build model
#' #==============
#' modelRandEff2 <- hmsc(formdata, niter = 200, nburn = 100, thin = 1, 
#' 						  verbose = FALSE)
#' 
#' #===============================
#' ### Construct correlation matrix
#' #===============================
#' corRes <- corRandomEff(modelRandEff2)
#' 
#' #=================================
#' ### Draw chord diagram with corRes
#' #=================================
#' library(circlize)
#' 
#' chordDiagram(corRes[, , 1, 1], symmetric = TRUE)
#' 
#' @keywords IO
#' @export
corRandomEff <-
function(hmsc,cor=TRUE,latentToUse=NULL,latentAutoToUse=NULL,burning=FALSE){
#### F. Guillaume Blanchet - April 2015, January 2016
##########################################################################################
	if(!any(names(hmsc$results$estimation)=="paramLatent" | names(hmsc$results$estimation)=="paramLatentAuto")){
		stop("Only works if 'paramLatent' was estimated (if there were random effects included in the  HMSC model)")
	}
	
	if(any(names(hmsc$results$estimation)=="paramLatent")){
		### Basic values
		nsp<-nrow(hmsc$results$estimation$paramLatent[[1,1]])
		nest<-nrow(hmsc$results$estimation$paramLatent)
		nRandom<-ncol(hmsc$results$estimation$paramLatent)
		
		### Check if latentToUse is too large
		if(!is.null(latentToUse)){
			if(burning){
				nburn<-nrow(hmsc$results$burning$paramLatent)
				nLatent<-matrix(NA,nrow=nest+nburn,ncol=nRandom)
				
				### Burning results
				for(i in 1:nRandom){
					for(j in 1:nburn){
						nLatent[j,i]<-ncol(hmsc$results$burning$paramLatent[[j,i]])
					}
				}
				
				### Estimation results
				for(i in 1:nRandom){
					for(j in 1:nest){
						nLatent[j+nburn,i]<-ncol(hmsc$results$estimation$paramLatent[[j,i]])
					}
				}
				
				minLatent<-min(nLatent)
				if(any(latentToUse>minLatent)){
					stop("The parameters referred to in 'latentToUse' are too large ")
				}
			}else{
				nLatent<-matrix(NA,nrow=nest,ncol=nRandom)
				
				for(i in 1:nRandom){
					for(j in 1:nest){
						nLatent[j,i]<-ncol(hmsc$results$estimation$paramLatent[[j,i]])
					}
				}
				
				minLatent<-min(nLatent)
				if(any(latentToUse>minLatent)){
					stop("The parameters referred to in 'latentToUse' are too large ")
				}
			}
		}
		
		### Include burning information
		if(burning){
			corMat<-array(dim=c(nsp,nsp,nest+nburn,nRandom))
			if(is.null(latentToUse)){
				for(i in 1:nRandom){
					for(j in 1:nburn){
						corMat[,,j,i]<-tcrossprod(hmsc$results$burn$paramLatent[[j,i]])
					}
					for(j in 1:nest){
						corMat[,,nburn+j,i]<-tcrossprod(hmsc$results$est$paramLatent[[j,i]])
					}
				}
			}else{
				for(i in 1:nRandom){
					for(j in 1:nburn){
						corMat[,,j,i]<-tcrossprod(hmsc$results$burn$paramLatent[[j,i]][,latentToUse])
					}
					for(j in 1:nest){
						corMat[,,nburn+j,i]<-tcrossprod(hmsc$results$est$paramLatent[[j,i]][,latentToUse])
					}
				}
			}
			dimnames(corMat)[[3]]<-c(rownames(hmsc$results$burn$paramLatent),rownames(hmsc$results$est$paramLatent))
			
		### Without burning information
		}else{
			corMat<-array(dim=c(nsp,nsp,nest,nRandom))
			if(is.null(latentToUse)){
				for(i in 1:nRandom){
					for(j in 1:nest){
						corMat[,,j,i]<-tcrossprod(hmsc$results$est$paramLatent[[j,i]])
					}
				}
			}else{
				for(i in 1:nRandom){
					for(j in 1:nest){
						corMat[,,j,i]<-tcrossprod(hmsc$results$est$paramLatent[[j,i]][,latentToUse])
					}
				}
			}
			dimnames(corMat)[[3]]<-rownames(hmsc$results$est$paramLatent)
		}
		
		### Convert to correlation
		if(cor==TRUE){
			niter<-dim(corMat)[3]
			for(i in 1:nRandom){
				for(j in 1:niter){
					corMat[,,j,i]<-cov2cor(corMat[,,j,i])
				}
			}
		}
		
		### Output
		dimnames(corMat)[[1]]<-dimnames(hmsc$results$estimation$paramLatent[[1,1]])[[1]]
		dimnames(corMat)[[2]]<-dimnames(hmsc$results$estimation$paramLatent[[1,1]])[[1]]
		if(nRandom==1){
			dimnames(corMat)[[4]]<-as.list(paste("randEff",1:nRandom,sep=""))
		}else{
			dimnames(corMat)[[4]]<-paste("randEff",1:nRandom,sep="")
		}
	}
	
	if(any(names(hmsc$results$estimation)=="paramLatentAuto")){
		### Basic values
		nsp<-nrow(hmsc$results$estimation$paramLatentAuto[[1,1]])
		nest<-nrow(hmsc$results$estimation$paramLatentAuto)
		nAuto<-ncol(hmsc$results$estimation$paramLatentAuto)
	
		### Check if latentAutoToUse is too large
		if(!is.null(latentAutoToUse)){
			if(burning){
				nburn<-nrow(hmsc$results$burning$paramLatentAuto)
				nLatentAuto<-matrix(NA,nrow=nest+nburn,ncol=nAuto)
				
				### Burning results
				for(i in 1:nAuto){
					for(j in 1:nburn){
						nLatentAuto[j,i]<-ncol(hmsc$results$burning$paramLatentAuto[[j,i]])
					}
				}
				
				### Estimation results
				for(i in 1:nAuto){
					for(j in 1:nest){
						nLatentAuto[j+nburn,i]<-ncol(hmsc$results$estimation$paramLatentAuto[[j,i]])
					}
				}
				
				minLatentAuto<-min(nLatentAuto)
				if(any(latentAutoToUse>minLatentAuto)){
					stop("The parameters referred to in 'latentAutoToUse' are too large ")
				}
			}else{
				nLatentAuto<-matrix(NA,nrow=nest,ncol=nAuto)
				
				for(i in 1:nAuto){
					for(j in 1:nest){
						nLatentAuto[j,i]<-ncol(hmsc$results$estimation$paramLatentAuto[[j,i]])
					}
				}
				
				minLatentAuto<-min(nLatentAuto)
				if(any(latentAutoToUse>minLatentAuto)){
					stop("The parameters referred to in 'latentAutoToUse' are too large ")
				}
			}
		}
		
		### Include burning information
		if(burning){
			corMatAuto<-array(dim=c(nsp,nsp,nest+nburn,nAuto))
			if(is.null(latentAutoToUse)){
				for(i in 1:nAuto){
					for(j in 1:nburn){
						corMatAuto[,,j,i]<-tcrossprod(hmsc$results$burn$paramLatentAuto[[j,i]])
					}
					for(j in 1:nest){
						corMatAuto[,,nburn+j,i]<-tcrossprod(hmsc$results$est$paramLatentAuto[[j,i]])
					}
				}
			}else{
				for(i in 1:nAuto){
					for(j in 1:nburn){
						corMatAuto[,,j,i]<-tcrossprod(hmsc$results$burn$paramLatentAuto[[j,i]][,latentAutoToUse])
					}
					for(j in 1:nest){
						corMatAuto[,,nburn+j,i]<-tcrossprod(hmsc$results$est$paramLatentAuto[[j,i]][,latentAutoToUse])
					}
				}
			}
			dimnames(corMatAuto)[[3]]<-c(rownames(hmsc$results$burn$paramLatentAuto),rownames(hmsc$results$est$paramLatentAuto))
			
		### Without burning information
		}else{
			corMatAuto<-array(dim=c(nsp,nsp,nest,nAuto))
			if(is.null(latentAutoToUse)){
				for(i in 1:nAuto){
					for(j in 1:nest){
						corMatAuto[,,j,i]<-tcrossprod(hmsc$results$est$paramLatentAuto[[j,i]])
					}
				}
			}else{
				for(i in 1:nAuto){
					for(j in 1:nest){
						corMatAuto[,,j,i]<-tcrossprod(hmsc$results$est$paramLatentAuto[[j,i]][,latentAutoToUse])
					}
				}
			}
			dimnames(corMatAuto)[[3]]<-rownames(hmsc$results$est$paramLatentAuto)
		}
		
		### Convert to correlation
		if(cor==TRUE){
			niter<-dim(corMatAuto)[3]
			for(i in 1:nAuto){
				for(j in 1:niter){
					corMatAuto[,,j,i]<-cov2cor(corMatAuto[,,j,i])
				}
			}
		}
		
		### Output
		dimnames(corMatAuto)[[1]]<-dimnames(hmsc$results$estimation$paramLatentAuto[[1,1]])[[1]]
		dimnames(corMatAuto)[[2]]<-dimnames(hmsc$results$estimation$paramLatentAuto[[1,1]])[[1]]
		if(nAuto==1){
			dimnames(corMatAuto)[[4]]<-as.list(paste("randEffAuto",1:nAuto,sep=""))
		}else{
			dimnames(corMatAuto)[[4]]<-paste("randEffAuto",1:nAuto,sep="")
		}
	}

	### Return results
	if(any(names(hmsc$results$estimation)=="paramLatent") & any(names(hmsc$results$estimation)=="paramLatentAuto")){
		res<-vector("list",length=2)
		res[[1]]<-corMat
		res[[2]]<-corMatAuto
		
		names(res)<-c("Random","Auto")
	}else{
		if(any(names(hmsc$results$estimation)=="paramLatent")){
			res<-corMat
		}else{
			if(any(names(hmsc$results$estimation)=="paramLatentAuto")){
				res<-corMatAuto
			}
		}
	}
	class(res)<-"corRandomEff"
	
	return(res)
}