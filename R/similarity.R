#' @title Community similarity induced by environmental similarity
#'
#' @description Calculates a community similarity resulting from environmental similarity.
#'
#' @param hmsc An object of the class \code{hmsc}.
#' @param gr A factor defining groups of samples. If \code{NULL}, all samples are considered as independent groups.
#'
#' @return
#'
#' A symmetrical covariance matrix between all pairs of samples
#'
#' @author F. Guillaume Blanchet
#'
#' @examples
#' #================
#' ### Generate data
#' #================
#' desc <- cbind(scale(1:50), scale(1:50)^2)
#' nspecies <- 20
#' commDesc <- communitySimul(X = desc, nsp = nspecies)
#'
#' #=============
#' ### Formatting
#' #=============
#' ### Format data
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE,
#' 						   interceptTr = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 2000, nburn = 1000, thin = 1, verbose = 100)
#'
#' #=======================
#' ### Calculate similarity
#' #=======================
#' similarity <- similarity(modelDesc)
#'
#' @keywords datagen
#' @keywords htest
#' @keywords univar
#' @keywords multivariate
#' @keywords regression
#' @export
similarity<-function(hmsc,gr=NULL){
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

	### Extract parameters
	param <- coef(hmsc)

	### Similarity if there are no traits or phylogeny
	if(!any(names(hmsc$results$estimation) %in% c("paramTr","paramPhylo"))){
		res <- X %*% cov(param$paramX) %*% t(X)
	}

	### Similarity if there are only traits
	if(any(names(hmsc$results$estimation) %in% "paramTr") & (!any(names(hmsc$results$estimation) %in% "paramPhylo"))){
		res <- X %*% (param$paramTr %*% cov(t(hmsc$data$Tr)) %*% t(param$paramTr)+param$varX) %*%t (X)
	}

	### Similarity if there is only phylogeny
	if(any(names(hmsc$results$estimation) %in% "paramPhylo") & (!any(names(hmsc$results$estimation) %in% "paramTrait"))){
		phylo <- hmsc$data$Phylo
		diag(phylo) <- NA
		phyloMean <- mean(phylo, na.rm = TRUE)

		res <- X %*% (tcrossprod(param$meansParamX)+(1-param$paramPhylo*phyloMean)*param$varX) %*%t (X)
	}

	### Similarity if there is only phylogeny
	if(any(names(hmsc$results$estimation) %in% c("paramTr","paramPhylo"))){
		phylo <- hmsc$data$Phylo
		diag(phylo) <- NA
		phyloMean <- mean(phylo, na.rm = TRUE)

		res <- X %*% (param$paramTr %*% cov(t(hmsc$data$Tr)) %*% t(param$paramTr)+(1-param$paramPhylo*phyloMean)*param$varX) %*%t (X)
	}

	### Return results
	return(res)
}
