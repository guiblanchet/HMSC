#' @title Coefficient of determination
#'
#' @description Calculate the coefficient of determination for a model
#'
#' @param hmsc An object of the class \code{hmsc}
#' @param newdata An optional object of class \code{HMSCdata} in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param adjust Logical. Whether an adjustement should be calculated on the calculation of the coefficient of determination. Default is \code{FALSE}.
#' @param averageSp Logical. Whether the coefficient of determination is calculated for all species independently (\code{FALSE}) or for the community as a whole (\code{TRUE}). Default is \code{TRUE}.
#'
#' @details
#'
#' If the community data is composed only of presence-absence data, this function calculates the coefficient of determination (\eqn{R^2}) proposed by Tjur (2009). The multivariate (or community-level) Tjur's \eqn{R^2} is calculated by averaging over the univariate (species-level) Tjur's \eqn{R^2}.
#'
#' If the HMSC model was built using a Gaussian model, the classical coefficient of determination (\eqn{R^2}) was calculated:
#'
#' \deqn{R^2 = 1 - \frac{SS_{resid}}{SS_{total}}}{1-SSresid/SStotal}.
#'
#' where \eqn{SS_{resid}}{SSresid} is the residual sum of squares and \eqn{SS_{total}}{SStotal} the total sum of squares.
#'
#' The adjustement to the coefficient of determination used here is the one proposed by Gelman and Pardoe (2006), which leads to a lower adjusted coefficient of determination than the classical adjustement because it accounts for uncertainty in the model variance. The adjustement used in this function is calculated as follows:
#'
#' \deqn{1 - \frac{n-3}{n-p-2}(1-R^2)}{1-((n-3)/n-p-2)(1-R2)}.
#'
#' @return
#'
#' If \code{averageSp} is \code{TRUE} a single value is returned leading to a community-level coefficient of determination.
#'
#' If \code{averageSp} is \code{FALSE} a vector with as many values as the number of species is returned leading to specie-level coefficient of determination.
#'
#' @references
#'
#' Tjur, T. (2009) Coefficients of determination in logistic regression models - A new proposal: The coefficient of discrimination. \emph{The American Statistician} \strong{63}, 366-372.
#'
#' Gelman, A. and I. Pardoe (2006) Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models. \emph{Technometrics} \strong{48}, 241-251.
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom stats predict
#' @examples
#'
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
#' formdata <- as.HMSCdata(Y = commDesc$data$Y, X = desc, interceptX = FALSE, interceptTr = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' modelDesc <- hmsc(formdata, niter = 2000, nburn = 1000, thin = 1, verbose = 100)
#'
#' #=======================================================
#' ### Calculate full coefficient of multiple determination
#' #=======================================================
#' Rsquared(modelDesc, averageSp = FALSE)
#'
#' @keywords univar, multivariate, regression
#' @export
Rsquared <- function(hmsc, newdata = NULL, adjust = FALSE, averageSp = TRUE) {
  ### Account for newdata
  if (is.null(newdata)) {
    Y <- hmsc$data$Y
  }
	if (!is.null(newdata)) {
	  Y <- newdata$Y
  }

  ### Warning for poisson and overPoisson models (this may change in the future)
  if (any(class(hmsc) == c("poisson", "overPoisson"))) {
      stop("R2 may only be calculated for probit and gaussian models... for now")
    }

  ### Number of species
  nsp <- ncol(Y)

  ### Calculate model estimates
  Ypred <- predict(hmsc, newdata=newdata)

  ### Probit model
  if (any(class(hmsc) == "probit")) {
    Y0 <- Y == 0
    Y1 <- Y == 1

    R2 <- numeric()

    for (i in 1:nsp) {
        R2[i] <- mean(Ypred[Y1[, i], i]) - mean(Ypred[Y0[, i], i])
    }
  }

  ### Gaussian model
  if (any(class(hmsc) == "gaussian")) {
    ### Total sums of squares per species
    ssY <- colSums(scale(Y,scale=FALSE)^2)

    ### Residual sums of squares per species
    ssRes <- colSums((Y-predict(hmsc))^2)

    ### Calculate R2
    R2 <- 1-ssRes/ssY
  }

  ### Adjustement
  if(adjust){
    nsite <- nrow(Y)

    #____________________________________________
    ### Count the number of explanatory variables
    #____________________________________________
    nexp <- 0
    ### X
    if(any(names(hmsc$data)=="X")){
      nexp <- nexp + ncol(hmsc$data$X)
    }

    ### Random
    if(any(names(hmsc$data)=="Random")){
      for(i in 1:ncol(hmsc$results$estimation$latent)){
        nexp <- nexp + max(sapply(hmsc$results$estimation$latent[,i],ncol))
      }
    }

    ### Auto
    if(any(names(hmsc$data)=="Auto")){
      for(i in 1:ncol(hmsc$results$estimation$latentAuto)){
        nexp <- nexp + max(sapply(hmsc$results$estimation$latentAuto[,i],ncol))
      }
    }

    ### Calculate adjusted R2
    R2 <- 1 - ((nsite-3)/(nsite - nexp - 2)) * (1 - R2)
  }

  ### Community-level R2
  if (averageSp) {
      R2 <- mean(R2)
  }

  return(R2)
}
