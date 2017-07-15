#' @title Convert MCMC parameters estimation to \code{\link{mcmc}} object
#'
#' @description Use the Markov Chain Monte Carlo results obtained from \code{\link{hmsc}} and convert them to an \code{\link{mcmc}} object.
#'
#' @param x An object of the class \code{hmsc}.
#' @param parameters A character string defining the parameters for which the convergions needs to be carried out.
#' @param burning Logical. Whether the burning iterations should be include (\code{TRUE}) or not (\code{FALSE}).
#' @param \dots Addtional arguments passed to \code{\link[coda]{mcmc}}.
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
#' @seealso \code{\link{mcmc}}, \code{\link{hmsc}}
#'
#' @importFrom coda mcmc
#' @importFrom coda as.mcmc
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
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE, interceptTr=FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 200, nburn = 100, thin = 1, verbose = FALSE)
#'
#' paramXmcmc <- as.mcmc(modelDesc, parameters = "paramX")
#'
#' @keywords IO
#' @export
as.mcmc.hmsc <-
function(x, parameters = "paramX", burning = FALSE, Random = 1, Auto = 1,...){
#### F. Guillaume Blanchet - April 2015, January 2016, June 2016
##########################################################################################

	### For latent
	if(parameters=="latent"){
		stop("an mcmc object for 'latent' should not be constructed")
	}

	### For paramLatent
	if(parameters == "paramLatent" | parameters=="paramLatentAuto"){
		nrandom <- ncol(x$results$estimation[[parameters]])
		nsp <- nrow(x$results$estimation[[parameters]][[1,1]])
		nlatent <- min(sapply(x$results$estimation[[parameters]],ncol))
		niter <- nrow(x$results$estimation$paramLatent)

		### Include burning information
		if(burning){
			nlatent<-min(sapply(x$results$estimation[[parameters]],ncol),sapply(x$results$burning[[parameters]],ncol))
			nburn <- length(x$results$burning[[parameters]])
			paramMCMC <- array(dim=c(nsp,nlatent,niter,nrandom))
			for(i in 1:nrandom){
				for(j in 1:nburn){
					paramMCMC[,,j,i] <- x$results$burning[[parameters]][[j,i]][,1:nlatent]
				}

				for(j in 1:niter){
					paramMCMC[,,j,i] <- x$results$estimation[[parameters]][[j,i]][,1:nlatent]
				}
			}

		### Without burning information
		}else{
			paramMCMC <- array(dim=c(nsp,nlatent,niter,nrandom))
			for(i in 1:nrandom){
				for(j in 1:niter){
					paramMCMC[,,j,i] <- x$results$estimation[[parameters]][[j,i]][,1:nlatent]
				}
			}
		}

		### Reorganize paramMCMC
		paramMCMCMat <- aperm(paramMCMC,c(3,1,2,4))
		dim(paramMCMCMat) <- c(niter,nsp*nlatent*nrandom)

		### Name the different dimensions of paramMCMCMat
		if(burning){
			rownames(paramMCMCMat) <- c(rownames(x$results$burning[[parameters]]),rownames(x$results$estimation[[parameters]]))
		}else{
			rownames(paramMCMCMat) <- rownames(x$results$estimation[[parameters]])
		}

		colNameRough<-expand.grid(colnames(x$data$Y), colnames(x$results$estimation[[parameters]][[1,1]])[1:nlatent], colnames(x$results$estimation[[parameters]]))

		if(nrow(colNameRough)>0){
			colnames(paramMCMCMat) <- paste(colNameRough[,1],".",colNameRough[,2],".",colNameRough[,3],sep="")
		}
	}

	### For varX
	if(parameters=="varX"){
		paramMCMC <- x$results$estimation[[parameters]]
		niter <- dim(paramMCMC)[3]

		lowerTri <- lower.tri(paramMCMC[,,1],diag=TRUE)
		lowerTriMatPointer <- which(lowerTri,arr.ind=TRUE)

		paramMCMCMat <- matrix(NA,niter,nrow(lowerTriMatPointer))
		for(i in 1:niter){
			paramMCMCMat[i,] <- paramMCMC[,,i][lowerTriMatPointer]
		}

		rownames(paramMCMCMat) <- dimnames(x$results$estimation$varX)[[3]]

		### Include burning information
		if(burning){
			paramBurnMCMC <- x$results$burning[[parameters]]
			nburn <- dim(paramBurnMCMC)[3]

			paramBurnMCMCMat <- matrix(NA,nburn,nrow(lowerTriMatPointer))
			for(i in 1:nburn){
				paramBurnMCMCMat[i,] <- paramMCMC[,,i][lowerTriMatPointer]
			}

			paramMCMCMat <- rbind(paramBurnMCMCMat,paramMCMCMat)
			rownames(paramMCMCMat) <- c(dimnames(x$results$burning$varX)[[3]],dimnames(x$results$estimation$varX)[[3]])
		}

		varXNames <- expand.grid(dimnames(x$results$estimation$varX)[[1]],dimnames(x$results$estimation$varX)[[1]])[which(lowerTri),]
		if(ncol(varXNames)>0){
			colnames(paramMCMCMat) <- paste(varXNames[,1],".",varXNames[,2],sep="")
		}
	}

	### meansParamX, varNormal, varPoisson, paramPhylo
	if(parameters=="meansParamX" | parameters=="varNormal" | parameters=="varPoisson" | parameters=="paramPhylo"){
		paramMCMCMat <- x$results$estimation[[parameters]]

		### Name rows and columns of matrix
		rownames(paramMCMCMat) <- rownames(x$results$estimation[[parameters]])
		colnames(paramMCMCMat) <- colnames(x$results$estimation[[parameters]])

		### Output
		if(burning){
			paramBurnMCMC <- x$results$burning[[parameters]]
			rownames(paramBurnMCMC) <- rownames(x$results$burning[[parameters]])
			colnames(paramBurnMCMC) <- colnames(x$results$burning[[parameters]])
			paramMCMCMat <- rbind(paramBurnMCMCMat,paramMCMCMat)
		}
	}

	### paramX
	if(parameters == "paramX" | parameters == "paramTr"){
		### Reorganize results
		paramMCMC <- x$results$estimation[[parameters]]
		paramMCMCMat <- aperm(paramMCMC,c(3,1,2))
		dim(paramMCMCMat) <- c(dim(x$results$estimation[[parameters]])[[3]],dim(x$results$estimation[[parameters]])[[1]]*dim(x$results$estimation[[parameters]])[[2]])

		### Name rows and columns of matrix
		rownames(paramMCMCMat) <- dimnames(x$results$estimation[[parameters]])[[3]]
		colNameRough <- expand.grid(dimnames(x$results$estimation[[parameters]])[[1]],dimnames(x$results$estimation[[parameters]])[[2]])

		if(nrow(colNameRough)>0){
			colnames(paramMCMCMat) <- paste(colNameRough[,1],".",colNameRough[,2],sep="")
		}

		### Include burning information
		if(burning){
			paramBurnMCMC <- x$results$burning[[parameters]]
			paramBurnMCMCMat <- aperm(paramBurnMCMC,c(3,1,2))
			dim(paramBurnMCMCMat) <- c(dim(x$results$burning[[parameters]])[[3]],dim(x$results$burning[[parameters]])[[1]]*dim(x$results$burning[[parameters]])[[2]])
			rownames(paramBurnMCMCMat) <- dimnames(x$results$burning[[parameters]])[[3]]
			paramMCMCMat <- rbind(paramBurnMCMCMat,paramMCMCMat)
		}
	}

	### Output
	res <- mcmc(paramMCMCMat, ...)

	return(res)
}
