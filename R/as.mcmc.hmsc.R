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
#' The parameters estimated for the latent variables (\code{paramLatent}) cannot be directly converted to an object of class \code{mcmc} because the procedure implemented in \code{\link{hmsc}} allows for the number of latent variables, and thus of parameters associated to these latent variables, to change. For this reason, if \code{paramLatent} is used, it is the matrix defining the correlations among species calculated from  \code{paramLatent} for which an \code{\link{mcmc}} object is constructed.
#'
#' The \code{\link{mcmc}} object associated to \code{paramLatent} includes information only for the correlations among single pairs of species, that is, only the lower triangle of the correlation matrix is considered.
#'
#' When there are multiple random effects estimated for the model, associated to each random effect there is a set of estimated latent parameters (\code{paramLatent}). For this reason, when the Markov Chain Monte Carlo results for \code{paramLatent} is converted with \code{as.mcmc}, parameters associated to each set of latent variables in \code{paramLatent} becomes an code{\link{mcmc}} object in a list. It is important to note that this "list" is \emph{NOT} an \code{mcmc.list}.
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
function(x,parameters="paramX",burning=FALSE,...){
#### F. Guillaume Blanchet - April 2015, January 2016, June 2016
##########################################################################################

	### For latent
	if(parameters=="latent"){
		stop("an mcmc object for 'latent' should not be constructed")
	}

	### For paramLatent
	if(parameters=="paramLatent" | parameters=="paramLatentAuto"){
		nrandom <- ncol(x$results$est[[parameters]])
		nsp <- nrow(x$results$est[[parameters]][[1,1]])
		niter <- nrow(x$results$est$paramLatent)

		### Include burning information
		if(burning){
			nburn <- length(x$results$burn[[parameters]])
			covMat <- array(dim=c(nsp,nsp,niter+nburn,nrandom))
			for(i in 1:nrandom){
				for(j in 1:nburn){
					covMat[,,j,i] <- tcrossprod(x$results$burn[[parameters]][[j,i]])
				}
				for(j in 1:niter){
					covMat[,,nburn+j,i] <- tcrossprod(x$results$est[[parameters]][[j,i]])
				}
			}

		### Without burning information
		}else{
			covMat <- array(dim=c(nsp,nsp,niter,nrandom))
			for(i in 1:nrandom){
				for(j in 1:niter){
					covMat[,,j,i] <- tcrossprod(x$results$est[[parameters]][[j,i]])
				}
			}
		}

		### Use only the lower triangle of the correlations matrices
		lowerTri <- lower.tri(covMat[,,1,1],diag=TRUE)
		lowerTriMatPointer <- which(lowerTri,arr.ind=TRUE)

		### Reorganize covMat with burning
		if(burning){
			paramMCMCMat <- array(dim=c(nburn+niter,nrow(lowerTriMatPointer),nrandom))
			for(i in 1:nrandom){
				for(j in 1:(nburn+niter)){
					paramMCMCMat[j,,i] <- covMat[,,j,i][lowerTriMatPointer]
				}
			}
			dimnames(paramMCMCMat)[[1]] <- c(names(x$results$burn[[parameters]]),names(x$results$est[[parameters]]))

		### Reorganize covMat without burning
		}else{
			paramMCMCMat <- array(dim=c(niter,nrow(lowerTriMatPointer),nrandom))
			for(i in 1:nrandom){
				for(j in 1:(niter)){
					paramMCMCMat[j,,i] <- covMat[,,j,i][lowerTriMatPointer]
				}
			}
			dimnames(paramMCMCMat)[[1]] <- names(x$results$est[[parameters]])
		}
		spNameRough <- expand.grid(rownames(x$results$est[[parameters]][[1,1]]),rownames(x$results$est[[parameters]][[1,1]]))[which(lowerTri),]
		dimnames(paramMCMCMat)[[2]] <- paste(spNameRough[,1],".",spNameRough[,2],sep="")
		dimnames(paramMCMCMat)[[3]] <- paste("randEff",1:dim(paramMCMCMat)[3])

		### Output
		res <- vector("list",length=nrandom)

		for(i in 1:nrandom){
			res[[i]] <- mcmc(paramMCMCMat[,,i], ...)
		}
	}else{
		### For varX
		if(parameters=="varX"){
			paramMCMC <- x$results$est[[parameters]]
			niter <- dim(paramMCMC)[3]

			lowerTri <- lower.tri(paramMCMC[,,1],diag=TRUE)
			lowerTriMatPointer <- which(lowerTri,arr.ind=TRUE)

			paramMCMCMat <- matrix(NA,niter,nrow(lowerTriMatPointer))
			for(i in 1:niter){
				paramMCMCMat[i,] <- paramMCMC[,,i][lowerTriMatPointer]
			}

			rownames(paramMCMCMat) <- dimnames(x$results$est$varX)[[3]]

			### Include burning information
			if(burning){
				paramBurnMCMC <- x$results$burn[[parameters]]
				nburn <- dim(paramBurnMCMC)[3]

				paramBurnMCMCMat <- matrix(NA,nburn,nrow(lowerTriMatPointer))
				for(i in 1:nburn){
					paramBurnMCMCMat[i,] <- paramMCMC[,,i][lowerTriMatPointer]
				}

				paramMCMCMat <- rbind(paramBurnMCMCMat,paramMCMCMat)
				rownames(paramMCMCMat) <- c(dimnames(x$results$burn$varX)[[3]],dimnames(x$results$est$varX)[[3]])
			}

			varXNames <- expand.grid(dimnames(x$results$est$varX)[[1]],dimnames(x$results$est$varX)[[1]])[which(lowerTri),]
			if(ncol(varXNames)>0){
				colnames(paramMCMCMat) <- paste(varXNames[,1],".",varXNames[,2],sep="")
			}

			### Output
			res <- mcmc(paramMCMCMat, ...)

		}else{
			if(parameters=="meansParamX" | parameters=="varNormal" | parameters=="paramPhylo"){
				paramMCMCMat <- x$results$est[[parameters]]

				### Name rows and columns of matrix
				rownames(paramMCMCMat) <- rownames(x$results$est[[parameters]])
				colnames(paramMCMCMat) <- colnames(x$results$est[[parameters]])

				### Output
				if(burning){
					paramBurnMCMC <- x$results$burn[[parameters]]
					rownames(paramBurnMCMC) <- rownames(x$results$burn[[parameters]])
					colnames(paramBurnMCMC) <- colnames(x$results$burn[[parameters]])
					paramMCMCMat <- rbind(paramBurnMCMCMat,paramMCMCMat)
				}
				res <- mcmc(paramMCMCMat, ...)

			}else{
				if(parameters=="paramPhylo"){

				}else{
					### Reorganize results
					paramMCMC <- x$results$est[[parameters]]
					paramMCMCMat <- aperm(paramMCMC,c(3,1,2))
					dim(paramMCMCMat) <- c(dim(x$results$est[[parameters]])[[3]],dim(x$results$est[[parameters]])[[1]]*dim(x$results$est[[parameters]])[[2]])

					### Name rows and columns of matrix
					rownames(paramMCMCMat) <- dimnames(x$results$est[[parameters]])[[3]]
					colNameRough <- expand.grid(dimnames(x$results$est[[parameters]])[[1]],dimnames(x$results$est[[parameters]])[[2]])

					if(nrow(colNameRough)>0){
						colnames(paramMCMCMat) <- paste(colNameRough[,1],".",colNameRough[,2],sep="")
					}

					### Include burning information
					if(burning){
						paramBurnMCMC <- x$results$burn[[parameters]]
						paramBurnMCMCMat <- aperm(paramBurnMCMC,c(3,1,2))
						dim(paramBurnMCMCMat) <- c(dim(x$results$burn[[parameters]])[[3]],dim(x$results$burn[[parameters]])[[1]]*dim(x$results$burn[[parameters]])[[2]])
						rownames(paramBurnMCMCMat) <- dimnames(x$results$burn[[parameters]])[[3]]
						paramMCMCMat <- rbind(paramBurnMCMCMat,paramMCMCMat)
					}

					### Output
					res <- mcmc(paramMCMCMat, ...)
				}
			}
		}
	}
	return(res)
}
