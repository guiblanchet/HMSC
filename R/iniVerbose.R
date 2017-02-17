#' @title Define what the function should print on screen
#'
#' @description Defines how the \code{verbose} arguments should be handled. This is meant to be an internal function.
#'
#' @param verbose a logical or numerical value
#' @param niter Number of of iterations carried out in the analysis
#'
#' @details
#'
#' If \code{verbose} is set to  \code{TRUE}, the number of iterations completed is printed on screen after every \code{niter}/5 iterations.
#' If \code{verbose} is a numerical value, a multiple of this value is printed on screen. 
#' If the number of iterations is smaller than 5, nothing is printed on screen.
#' 
#' @author F. Guillaume Blanchet
#' @keywords print
#' @export
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
