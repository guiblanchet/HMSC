iniVerbose <-
function(verbose,niter){
	if(is.logical(verbose)){
		if(verbose){
			### If verbose = TRUE, the number of iterations completed is printed on screen every 500 iterations
			verbose<-niter/5
		}else{
			### This ensure that nothing is printed on screen
			verbose<-niter+5
		}
	}
	### When verbose is a numerical value
	if(niter<5){
		verbose<-niter+5
	}else{
		verbose<-round(verbose)
	}
	return(verbose)
}
