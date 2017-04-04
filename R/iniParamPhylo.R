#' @importFrom stats runif
iniParamPhylo <-
function(priors,paramPhylo=NULL){

	if(is.null(paramPhylo)){
		paramPhylo<-sample(priors$param$paramPhylo[,1],1)
	}else{
		paramPhylo<-priors$param$paramPhylo[which.min(abs(abs(priors$param$paramPhylo[,1])-paramPhylo)),1]
	}

	### List of priors
	param<-list(paramPhylo=paramPhylo)

	names(param$paramPhylo)<-"phylo"

	return(param)
}
