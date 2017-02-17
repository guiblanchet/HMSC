#' @rdname nameResult
#' @export
nameXTr<-function(data,result,niter,nburn,thin){
	if(length(result)==2){
		#----------
		### Burning
		#----------
		### paramX
		dimnames(result[[1]][[1]])[[1]]<-colnames(data$Y)
		dimnames(result[[1]][[1]])[[3]]<-paste("iter",seq(1,nburn,by=thin),sep="")
		dimnames(result[[1]][[1]])[[2]]<-colnames(data$X)

		### paramTr
		dimnames(result[[1]][[2]])[[3]]<-paste("iter",seq(1,nburn,by=thin),sep="")
		dimnames(result[[1]][[2]])[[1]]<-colnames(data$X)
		dimnames(result[[1]][[2]])[[2]]<-rownames(data$Tr)
		
		### varX
		dimnames(result[[1]][[3]])[[3]]<-paste("iter",seq(1,nburn,by=thin),sep="")
		dimnames(result[[1]][[3]])[[1]]<-colnames(data$X)
		dimnames(result[[1]][[3]])[[2]]<-colnames(data$X)

		#-------------
		### Estimation
		#-------------
		### paramX
		dimnames(result[[2]][[1]])[[1]]<-colnames(data$Y)
		dimnames(result[[2]][[1]])[[3]]<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		dimnames(result[[2]][[1]])[[2]]<-colnames(data$X)

		### paramTr
		dimnames(result[[2]][[2]])[[3]]<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		dimnames(result[[2]][[2]])[[1]]<-colnames(data$X)
		dimnames(result[[2]][[2]])[[2]]<-rownames(data$Tr)
		
		### varX
		dimnames(result[[2]][[3]])[[3]]<-paste("iter",seq(nburn+1,niter,by=thin),sep="")
		dimnames(result[[2]][[3]])[[1]]<-colnames(data$X)
		dimnames(result[[2]][[3]])[[2]]<-colnames(data$X)
	}
	
	if(length(result)==1){
		### paramX
		dimnames(result[[1]][[1]])[[1]]<-colnames(data$Y)
		dimnames(result[[1]][[1]])[[3]]<-paste("iter",seq(1,niter,by=thin),sep="")
		dimnames(result[[1]][[1]])[[2]]<-colnames(data$X)

		### paramTr
		dimnames(result[[1]][[2]])[[1]]<-paste("iter",seq(1,niter,by=thin),sep="")
		dimnames(result[[1]][[2]])[[1]]<-colnames(data$X)
		dimnames(result[[1]][[2]])[[2]]<-rownames(data$Tr)
		
		### varX
		dimnames(result[[1]][[3]])[[3]]<-paste("iter",seq(1,niter,by=thin),sep="")
		dimnames(result[[1]][[3]])[[1]]<-colnames(data$X)
		dimnames(result[[1]][[3]])[[2]]<-colnames(data$X)
	}

	return(result)
}
