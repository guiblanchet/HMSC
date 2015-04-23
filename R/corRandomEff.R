corRandomEff <-
function(hmsc,burning=FALSE){
#### F. Guillaume Blanchet - April 2015
##########################################################################################
	if(!any(names(hmsc$est)=="paramLatent")){
		stop("Only works if 'paramLatent' was estimated (if there were random effects included in the  HMSC model)")
	}
	
	### Basic values
	nsp<-nrow(hmsc$burn$paramX)
	niter<-length(hmsc$est$paramLatent)
	nrandom<-length(hmsc$est$paramLatent[[1]])
	### Include burning information
	if(burning){
		nburn<-length(hmsc$burn$paramLatent)
		corMat<-array(dim=c(nsp,nsp,niter+nburn,nrandom))
		for(i in 1:nrandom){
			for(j in 1:nburn){
				corMat[,,j,i]<-cov2cor(tcrossprod(hmsc$burn$paramLatent[[j]][[i]])+diag(nsp))
			}
			for(j in 1:niter){
				corMat[,,nburn+j,i]<-cov2cor(tcrossprod(hmsc$est$paramLatent[[j]][[i]])+diag(nsp))
			}
		}
	
		dimnames(corMat)[[3]]<-c(names(hmsc$burn$paramLatent),names(hmsc$est$paramLatent))

	### Without burning information
	}else{
		corMat<-array(dim=c(nsp,nsp,niter,nrandom))
		for(i in 1:nrandom){
			for(j in 1:niter){
					corMat[,,j,i]<-cov2cor(tcrossprod(hmsc$est$paramLatent[[j]][[i]])+diag(nsp))
			}
		}

		dimnames(corMat)[[3]]<-names(hmsc$est$paramLatent)
	}
	
	### Output
	dimnames(corMat)[[1]]<-dimnames(hmsc$est$paramX)[[1]]
	dimnames(corMat)[[2]]<-dimnames(hmsc$est$paramX)[[1]]
	dimnames(corMat)[[4]]<-paste("randEff",1:nrandom,sep="")
	
	class(corMat)<-"corRandomEff"
	
	return(corMat)
}