as.mcmc <-
function(hmsc,parameters="paramX",burning=FALSE){
#### F. Guillaume Blanchet - April 2015
##########################################################################################
	
	### For latent
	if(parameters=="latent"){
		stop("an mcmc object for 'latent' should not be constructed")
	}
	
	### For paramLatent
	if(parameters=="paramLatent"){
		nrandom<-length(hmsc$est[[parameters]][[1]])
		nsp<-nrow(hmsc$est[[parameters]][[1]][[1]])
		niter<-length(hmsc$est$paramLatent)
		
		### Include burning information
		if(burning){
			nburn<-length(hmsc$burn$paramLatent)
			corMat<-array(dim=c(nsp,nsp,niter+nburn,nrandom))
			for(i in 1:nrandom){
				for(j in 1:nburn){
					corMat[,,j,i]<-cov2cor(tcrossprod(hmsc$burn[[parameters]][[j]][[i]])+diag(nsp))
				}
				for(j in 1:niter){
					corMat[,,nburn+j,i]<-cov2cor(tcrossprod(hmsc$est[[parameters]][[j]][[i]])+diag(nsp))
				}
			}
		
		### Without burning information
		}else{
			corMat<-array(dim=c(nsp,nsp,niter,nrandom))
			for(i in 1:nrandom){
				for(j in 1:niter){
					corMat[,,j,i]<-cov2cor(tcrossprod(hmsc$est[[parameters]][[j]][[i]])+diag(nsp))
				}
			}
		}
		
		### Use only the lower triangle of the correlations matrices
		lowerTri<-row(corMat[,,1,1]) > col(corMat[,,1,1])
		lowerTriMatPointer<-which(lowerTri,arr.ind=TRUE)
		
		### Reorganize corMat with burning
		if(burning){
			paramMCMCMat<-array(dim=c(nburn+niter,nrow(lowerTriMatPointer),nrandom))
			for(i in 1:nrandom){
				for(j in 1:(nburn+niter)){
					paramMCMCMat[j,,i]<-corMat[,,j,i][lowerTriMatPointer]
				}
			}
			dimnames(paramMCMCMat)[[1]]<-c(names(hmsc$burn$paramLatent),names(hmsc$est$paramLatent))
		
		### Reorganize corMat without burning
		}else{
			paramMCMCMat<-array(dim=c(niter,nrow(lowerTriMatPointer),nrandom))
			for(i in 1:nrandom){
				for(j in 1:(niter)){
					paramMCMCMat[j,,i]<-corMat[,,j,i][lowerTriMatPointer]
				}
			}
			dimnames(paramMCMCMat)[[1]]<-names(hmsc$est$paramLatent)
		}
		spNameRough<-expand.grid(dimnames(hmsc$est$paramX)[[1]],dimnames(hmsc$est$paramX)[[1]])[which(lowerTri),]
		dimnames(paramMCMCMat)[[2]]<-paste(spNameRough[,1],".",spNameRough[,2],sep="")
		dimnames(paramMCMCMat)[[3]]<-paste("randEff",1:dim(paramMCMCMat)[3])
		
		### Output
		res<-vector("list",length=nrandom)
		
		for(i in 1:nrandom){
			res[[i]]<-mcmc(paramMCMCMat[,,i])
		}
	}else{
		### For varX
		if(parameters=="varX"){
			paramMCMC<-hmsc$est[[parameters]]
			niter<-dim(paramMCMC)[3]
			
			lowerTri<-row(paramMCMC[,,1]) > col(paramMCMC[,,1])
			lowerTriMatPointer<-which(lowerTri,arr.ind=TRUE)
			
			paramMCMCMat<-matrix(NA,niter,nrow(lowerTriMatPointer))
			for(i in 1:niter){
				paramMCMCMat[i,]<-paramMCMC[,,i][lowerTriMatPointer]
			}
			
			rownames(paramMCMCMat)<-dimnames(hmsc$est$varX)[[3]]
			
			### Include burning information
			if(burning){
				paramBurnMCMC<-hmsc$burn[[parameters]]
				nburn<-dim(paramBurnMCMC)[3]
				
				paramBurnMCMCMat<-matrix(NA,nburn,nrow(lowerTriMatPointer))
				for(i in 1:nburn){
					paramBurnMCMCMat[i,]<-paramMCMC[,,i][lowerTriMatPointer]
				}
				
				paramMCMCMat<-rbind(paramBurnMCMCMat,paramMCMCMat)
				rownames(paramMCMCMat)<-c(dimnames(hmsc$burn$varX)[[3]],dimnames(hmsc$est$varX)[[3]])
			}
			
			varXNames<-expand.grid(dimnames(hmsc$est$varX)[[1]],dimnames(hmsc$est$varX)[[1]])[which(lowerTri),]
			colnames(paramMCMCMat)<-paste(varXNames[,1],".",varXNames[,2],sep="")
			
			### Output
			res<-mcmc(paramMCMCMat)
			
		}else{
			### Reorganize results
			paramMCMC<-hmsc$est[[parameters]]
			paramMCMCMat<-aperm(paramMCMC,c(3,1,2))
			dim(paramMCMCMat)<-c(dim(hmsc$est[[parameters]])[[3]],dim(hmsc$est[[parameters]])[[1]]*dim(hmsc$est[[parameters]])[[2]])
		
			### Name rows and columns of matrix
			rownames(paramMCMCMat)<-dimnames(hmsc$est[[parameters]])[[3]]
			colNameRough<-expand.grid(dimnames(hmsc$est[[parameters]])[[1]],dimnames(hmsc$est[[parameters]])[[2]])
			colnames(paramMCMCMat)<-paste(colNameRough[,1],".",colNameRough[,2],sep="")
			
			### Include burning information
			if(burning){
				paramBurnMCMC<-hmsc$burn[[parameters]]
				paramBurnMCMCMat<-aperm(paramBurnMCMC,c(3,1,2))
				dim(paramBurnMCMCMat)<-c(dim(hmsc$burn[[parameters]])[[3]],dim(hmsc$burn[[parameters]])[[1]]*dim(hmsc$burn[[parameters]])[[2]])
				rownames(paramBurnMCMCMat)<-dimnames(hmsc$burn[[parameters]])[[3]]
				paramMCMCMat<-rbind(paramBurnMCMCMat,paramMCMCMat)
			}
			
			### Output
			res<-mcmc(paramMCMCMat)
		}
	}
	
	return(res)
}