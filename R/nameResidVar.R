nameResidVar<-function(data,result,niter,nburn,thin,listLev=4){
### listLev = part of the list to look into
	if(length(result)==2){
		#----------
		### Burning
		#----------
		seqBurn<-seq(1,nburn,by=thin)

		rownames(result[[1]][[listLev]])<-paste("iter",seqBurn,sep="")
		colnames(result[[1]][[listLev]])<-colnames(data$Y)

		#-------------
		### Estimation
		#-------------
		seqIter<-seq(nburn+1,niter,by=thin)

		rownames(result[[2]][[listLev]])<-paste("iter",seqIter,sep="")
		colnames(result[[2]][[listLev]])<-colnames(data$Y)
	}

	if(length(result)==1){
		seqIter<-seq(1,nburn,by=thin)

		rownames(result[[1]][[listLev]])<-paste("iter",seqIter,sep="")
		colnames(result[[1]][[listLev]])<-colnames(data$Y)
	}

	return(result)
}
