#' @title Coefficient of determination
#'
#' @description Calculate the coefficient of determination for a model
#'
#' @param hmsc An object of the class \code{hmsc}
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
Rsquared <- function(hmsc, averageSp = TRUE) {
    ### Number of species
    nsp <- ncol(hmsc$data$Y)

    ### Calculate model estimates
    Ypred <- predict(hmsc)

    ### Probit model
    if (any(class(hmsc) == "probit")) {
      Y0 <- hmsc$data$Y == 0
      Y1 <- hmsc$data$Y == 1

      R2 <- numeric()

      for (i in 1:nsp) {
          R2[i] <- mean(Ypred[Y1[, i], i]) - mean(Ypred[Y0[, i], i])
      }
    }

    ### Gaussian model
    if (any(class(hmsc) == "gaussian")) {
      ### Total sums of squares per species
      ssY <- colSums(scale(hmsc$data$Y,scale=FALSE)^2)

      ### Residual sums of squares per species
      ssRes <- colSums((hmsc$data$Y-predict(hmsc))^2)

      ### Calculate R2
      R2 <- 1-ssRes/ssY
    }

    ### Community-level R2
    if (averageSp) {
        R2 <- mean(R2)
    }

    return(R2)
}
