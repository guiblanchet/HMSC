nameAuto<-function(data,result,niter,nburn,thin,listLev=1){
### listLev = part of the list to look into
	if(length(result)==2){
		#----------
		### Burning
		#----------
		seqIter<-seq(1,nburn,by=thin)

		### paramLatentAuto
		rownames(result[[1]][[listLev]])<-paste("iter",seq(1,nburn,by=thin),sep="")
		colnames(result[[1]][[listLev]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[1]][[listLev]])){
			for(j in 1:ncol(result[[1]][[listLev]])){
				rownames(result[[1]][[listLev]][[i,j]])<-colnames(data$Y)
				colnames(result[[1]][[listLev]][[i,j]])<-paste("latentAuto",1:ncol(result[[1]][[listLev]][[i,j]]),sep="")
			}
		}

		### latentAuto
		rownames(result[[1]][[listLev+1]])<-paste("iter",seq(1,nburn,by=thin),sep="")
		colnames(result[[1]][[listLev+1]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[1]][[listLev+1]])){
			for(j in 1:ncol(result[[1]][[listLev+1]])){
				rownames(result[[1]][[listLev+1]][[i,j]])<-levels(data$Auto[[j]][,1])
				colnames(result[[1]][[listLev+1]][[i,j]])<-paste("latentAuto",1:ncol(result[[1]][[listLev+1]][[i,j]]),sep="")
			}
		}

		### paramAuto
		rownames(result[[1]][[listLev+2]])<-paste("iter",seq(1,nburn,by=thin),sep="")
		colnames(result[[1]][[listLev+2]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[1]][[listLev+2]])){
			for(j in 1:ncol(result[[1]][[listLev+2]])){
				rownames(result[[1]][[listLev+2]][[i,j]])<-paste("latentAuto",1:nrow(result[[1]][[listLev+2]][[i,j]]),sep="")
				colnames(result[[1]][[listLev+2]][[i,j]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")
			}
		}

		#-------------
		### Estimation
		#-------------
		### paramLatentAuto
		rownames(result[[2]][[listLev]])<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		colnames(result[[2]][[listLev]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[2]][[listLev]])){
			for(j in 1:ncol(result[[2]][[listLev]])){
				rownames(result[[2]][[listLev]][[i,j]])<-colnames(data$Y)
				colnames(result[[2]][[listLev]][[i,j]])<-paste("latentAuto",1:ncol(result[[2]][[listLev]][[i,j]]),sep="")
			}
		}
		### latentAuto
		rownames(result[[2]][[listLev+1]])<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		colnames(result[[2]][[listLev+1]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[2]][[listLev+1]])){
			for(j in 1:ncol(result[[2]][[listLev+1]])){
				rownames(result[[2]][[listLev+1]][[i,j]])<-levels(data$Auto[[j]][,1])
				colnames(result[[2]][[listLev+1]][[i,j]])<-paste("latentAuto",1:ncol(result[[2]][[listLev+1]][[i,j]]),sep="")
			}
		}

		### paramAuto
		rownames(result[[2]][[listLev+2]])<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		colnames(result[[2]][[listLev+2]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[2]][[listLev+2]])){
			for(j in 1:ncol(result[[2]][[listLev+2]])){
				rownames(result[[2]][[listLev+2]][[i,j]])<-paste("latentAuto",1:nrow(result[[2]][[listLev+2]][[i,j]]),sep="")
				colnames(result[[2]][[listLev+2]][[i,j]])<-paste("Auto",1:ncol(result[[2]][[listLev]]),sep="")
			}
		}
	}

	if(length(result)==1){
		### paramLatentAuto
		rownames(result[[1]][[listLev]])<-paste("iter",seq(1,niter,by=thin),sep="")
		colnames(result[[1]][[listLev]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[1]][[listLev]])){
			for(j in 1:ncol(result[[1]][[listLev]])){
				rownames(result[[1]][[listLev]][[i,j]])<-colnames(data$Y)
				colnames(result[[1]][[listLev]][[i,j]])<-paste("latentAuto",1:ncol(result[[1]][[listLev]][[i,j]]),sep="")
			}
		}

		### latentAuto
		rownames(result[[1]][[listLev+1]])<-paste("iter",seq(1,niter,by=thin),sep="")
		colnames(result[[1]][[listLev+1]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[1]][[listLev+1]])){
			for(j in 1:ncol(result[[1]][[listLev+1]])){
				rownames(result[[1]][[listLev+1]][[i,j]])<-levels(data$Auto[[j]][,1])
				colnames(result[[1]][[listLev+1]][[i,j]])<-paste("latentAuto",1:ncol(result[[1]][[listLev+1]][[i,j]]),sep="")
			}
		}

		### paramAuto
		rownames(result[[1]][[listLev+2]])<-paste("iter",seq(1,niter,by=thin),sep="")
		colnames(result[[1]][[listLev+2]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")

		for(i in 1:nrow(result[[1]][[listLev+2]])){
			for(j in 1:ncol(result[[1]][[listLev+2]])){
				rownames(result[[1]][[listLev+2]][[i,j]])<-paste("latentAuto",1:nrow(result[[1]][[listLev+2]][[i,j]]),sep="")
				colnames(result[[1]][[listLev+2]][[i,j]])<-paste("Auto",1:ncol(result[[1]][[listLev]]),sep="")
			}
		}
	}

	return(result)
}
