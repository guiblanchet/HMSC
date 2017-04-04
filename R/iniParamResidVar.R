iniParamResidVar <-
function(data,residVar=NULL){

	### Number of species
	nsp<-ncol(data$Y)
	if(is.null(residVar)){
		residVar<-rep(1,nsp)
	}else{
		if(length(residVar)!=nsp){
			stop("'residVar' is defined for a number of species that is different the one in 'data$Y'")
		}
	}

	### List of priors
	param<-list(residVar=residVar)

	return(param)
}
