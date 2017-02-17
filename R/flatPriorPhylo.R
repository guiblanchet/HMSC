#' @rdname flatPriorX
#' @export
flatPriorPhylo <-
function(paramPhylo=NULL){

	if(is.null(paramPhylo)){
		paramPhylo<-cbind(round(seq(0,1,by=0.01),2),1/101) # 101 is the number of entries in the grid
		colnames(paramPhylo)<-c("grid","weight")
	}
	
	### List of priors
	priors<-paramPhylo
	
	return(priors)
}
