#' @title Community similarity induced by environmental similarity
#'
#' @description Calculates a community similarity resulting from environmental similarity.
#'
#' @param hmsc An object of the class \code{hmsc}.
#' @param gr A factor defining groups of samples. If \code{NULL}, all samples are considered as independent groups.
#'
#' @details
#'
#' @return
#'
#' A symmetrical covariance matrix
#'
#' @author F. Guillaume Blanchet
#'
#' @examples
#'
#' @keywords datagen
#' @keywords htest
#' @keywords univar
#' @keywords multivariate
#' @keywords regression
#' @export
hmscSimil<-function(hmsc,gr=NULL){
#### F. Guillaume Blanchet
##########################################################################################
	### General Check
	if(!any(names(hmsc$results[[1]]) == "paramX")){
		stop("The HMSC model does not include any 'X'")
	}

	if(!is.factor(gr)){
		stop("'gr' needs to be a factor")
	}

	### Group the sites
	if(!is.null(gr)){
		nGr <- nlevels(gr)
		X <- matrix(NA, nrow = nGr, ncol = ncol(hmsc$data$X))

		for (i in 1:nGr) {
  		X[,i] <- tapply(hmsc$data$X[,i], gr, mean)
		}
	}else{
		X <- hmsc$data$X
	}

	### Extract paramX
	paramX <- modelDesc$results$estimation$paramX

	### Similarity if there are no traits or phylogeny
	if(!any(names(hmsc$results[[1]]) %in% c("paramTr","paramPhylo"))){
		res <- X[1,] %*% cov(paramX[,,1])%*%t(X[2,])
	}
 ### ICI ###


}
