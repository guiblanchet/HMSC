#' @rdname nameResult
#' @export
nameRandom<-function(data,result,niter,nburn,thin,listLev=1){
### listLev = part of the list to look into
	if(length(result)==2){
		#----------
		### Burning
		#----------
		seqIter<-seq(1,nburn,by=thin)

		### paramLatent
		rownames(result[[1]][[listLev]])<-paste("iter",seq(1,nburn,by=thin),sep="")
		colnames(result[[1]][[listLev]])<-paste("Random",1:ncol(result[[1]][[listLev]]),sep="")
		
		for(i in 1:nrow(result[[1]][[listLev]])){
			for(j in 1:ncol(result[[1]][[listLev]])){
				rownames(result[[1]][[listLev]][[i,j]])<-colnames(data$Y)
				colnames(result[[1]][[listLev]][[i,j]])<-paste("latent",1:ncol(result[[1]][[listLev]][[i,j]]),sep="")
			}
		}
		
		### latent
		rownames(result[[1]][[listLev+1]])<-paste("iter",seq(1,nburn,by=thin),sep="")
		colnames(result[[1]][[listLev+1]])<-paste("Random",1:ncol(result[[1]][[listLev]]),sep="")
		
		for(i in 1:nrow(result[[1]][[listLev+1]])){
			for(j in 1:ncol(result[[1]][[listLev+1]])){
				rownames(result[[1]][[listLev+1]][[i,j]])<-levels(data$Random[,j])
				colnames(result[[1]][[listLev+1]][[i,j]])<-paste("latent",1:ncol(result[[1]][[listLev+1]][[i,j]]),sep="")
			}
		}

		#-------------
		### Estimation
		#-------------
		### paramLatent
		rownames(result[[2]][[listLev]])<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		colnames(result[[2]][[listLev]])<-paste("Random",1:ncol(result[[1]][[listLev]]),sep="")
		
		for(i in 1:nrow(result[[2]][[listLev]])){
			for(j in 1:ncol(result[[2]][[listLev]])){
				rownames(result[[2]][[listLev]][[i,j]])<-colnames(data$Y)
				colnames(result[[2]][[listLev]][[i,j]])<-paste("latent",1:ncol(result[[2]][[listLev]][[i,j]]),sep="")
			}
		}
		### latent
		rownames(result[[2]][[listLev+1]])<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		colnames(result[[2]][[listLev+1]])<-paste("Random",1:ncol(result[[1]][[listLev]]),sep="")
		
		for(i in 1:nrow(result[[2]][[listLev+1]])){
			for(j in 1:ncol(result[[2]][[listLev+1]])){
				rownames(result[[2]][[listLev+1]][[i,j]])<-levels(data$Random[,j])
				colnames(result[[2]][[listLev+1]][[i,j]])<-paste("latent",1:ncol(result[[2]][[listLev+1]][[i,j]]),sep="")
			}
		}
	}

	if(length(result)==1){
		### paramLatent
		rownames(result[[1]][[listLev]])<-paste("iter",seq(1,niter,by=thin),sep="")
		colnames(result[[1]][[listLev]])<-paste("Random",1:ncol(result[[1]][[listLev]]),sep="")
		
		for(i in 1:nrow(result[[1]][[listLev]])){
			for(j in 1:ncol(result[[1]][[listLev]])){
				rownames(result[[1]][[listLev]][[i,j]])<-colnames(data$Y)
				colnames(result[[1]][[listLev]][[i,j]])<-paste("latent",1:ncol(result[[1]][[listLev]][[i,j]]),sep="")
			}
		}
		
		### latent
		rownames(result[[1]][[listLev+1]])<-paste("iter",seq(1,niter,by=thin),sep="")
		colnames(result[[1]][[listLev+1]])<-paste("Random",1:ncol(result[[1]][[listLev]]),sep="")
		
		for(i in 1:nrow(result[[1]][[listLev+1]])){
			for(j in 1:ncol(result[[1]][[listLev+1]])){
				rownames(result[[1]][[listLev+1]][[i,j]])<-nlevels(data$Random[,j])
				colnames(result[[1]][[listLev+1]][[i,j]])<-paste("latent",1:ncol(result[[1]][[listLev+1]][[i,j]]),sep="")
			}
		}
	}
	
	return(result)
}
