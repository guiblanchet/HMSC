#' @title Convert lists of hmsc objects to \code{\link{mcmc.list}} object
#'
#' @description Use the Markov Chain Monte Carlo results obtained from multiple \code{\link{hmsc}} and convert them to an \code{\link{mcmc.list}} object.
#'
#' @param x A list where each level contains an object of the class \code{hmsc}.
#' @param parameters A character string defining the parameters for which the \code{\link{mcmc.list}} needs to be constructed.
#' @param burning Logical. Whether the burning iterations should be include (\code{TRUE}) or not (\code{FALSE}).
#' @param \dots Additional arguments passed to \code{\link[coda]{mcmc.list}}.
#'
#' @details
#'
#' Since the algorithm used adapts the number of latent variables to use, the \code{\link{mcmc}} object associated to \code{paramLatent} only includes information associated to the latent variables that are available for all iterations.
#'
#' @return
#'
#' An object of class \code{\link{mcmc}}.
#'
#' @author F. Guillaume Blanchet
#'
#' @seealso \code{\link{mcmc.list}}, \code{\link{hmsc}}
#'
#' @importFrom coda mcmc.list
#' @examples
#'
#' #================
#' ### Generate data
#' #================
#' desc <- cbind(1, scale(1:50), scale(1:50)^2)
#' nspecies <- 20
#' commDesc <- communitySimul(X = desc, nsp = nspecies)
#'
#' #=============
#' ### Formatting
#' #=============
#' ### Format data
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#'
#' model <- vector("list",length=5)
#'
#' for(i in 1:length(model)){
#' 	model[[i]] <- hmsc(formdata, niter = 200, nburn = 100, thin = 1, verbose = FALSE)
#' }
#' paramXmcmc <- as.mcmc.list(model, parameters = "paramX")
#'
#' @keywords IO
#' @method as.mcmc.list hmsc
#' @export
as.mcmc.list.hmsc <-
function(x, parameters = NULL, burning = FALSE, ...){
#### F. Guillaume Blanchet - July 2017
##########################################################################################
	if(is.null(parameters)){
		stop("'parameters' must be specified by the user")
	}

	### Checks
	if(!is.list(x)){
		stop("'x' needs to be a list")
	}

	xClass <- sapply(x, function(x) class(x)[1])
	if(!all(xClass=="hmsc")){
		stop("all levels of 'x' needs to be of class 'hmsc'")
	}

	nrep <- length(x)
	paramMCMC <- vector("list",length=nrep)
	for(i in 1:nrep){
		paramMCMC[[i]] <- as.mcmc(x[[i]], parameters = parameters, burning = burning)
	}

	### Output
	res <- as.mcmc.list(paramMCMC)
	return(res)
}
