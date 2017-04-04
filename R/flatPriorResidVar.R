flatPriorResidVar <- 
function(residVar=NULL){
	### Prior for residVar
	if(is.null(residVar)){
		residVar<-c(1,0.3) # shape and scale
	}

	### List of priors
	priors<-list(residVar=residVar)

	return(priors)
}
