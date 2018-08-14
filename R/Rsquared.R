#' @title Coefficient of determination
#'
#' @description Calculate the coefficient of determination for a model
#'
#' @param hmsc An object of the class \code{hmsc}
#' @param newdata An optional object of class \code{HMSCdata} in which to look for variables with which to predict. If omitted, the fitted values are used.
#' @param type A character string describing the type of coefficient of determination to calculate. (See details).
#' @param adjust Logical. Whether an adjustement should be calculated on the calculation of the coefficient of determination. Default is \code{FALSE}.
#' @param averageSp Logical. Whether the coefficient of determination is calculated for all species independently (\code{FALSE}) or for the community as a whole (\code{TRUE}). Default is \code{TRUE}.
#' @param indepentSite Logical. Whether the coefficient of determination is calculated for each site independently (\code{TRUE}) or for all sites together (\code{FALSE}). This argument is only active when McFadden's pseudo-\eqn{R^2} is used. Default is \code{FALSE}.
#'
#' @details
#'
#' Four different version (\code{type}) are currently available in this function, "tjur" (Tjur 2009) which can only be used for binary (0, 1) data and "ols" which can only be used for Gaussian model. As for "efron" (Efron 1978) and "nakagawa" (Nakagawa and Schielzeth 2013), they can be used for all types of models.
#'
#' The multivariate (or community-level) \eqn{R^2}{R2} is calculated by averaging over the univariate (species-level) \eqn{R^2}{R2}.
#'
#' The \code{type = "ols"} can only be used with Gaussian models, the classical coefficient of determination (\eqn{R^2}{R2}) is calculated:
#'
#' \deqn{R^2 = 1 - \frac{(y_i-\hat{y})_i^2}{(y_i-\bar{y})^2}}{R2 = 1-((y-yhat)^2)/((y-ybar)^2)}.
#'
#' where \eqn{y_i}{y} is a response variable (species) at site i, \eqn{\hat{y}_i}{yhat} is the model predicted value for site i and \eqn{\bar{y}_i}{ybar} is the average of the response variable (species).
#'
#' The \code{type = "tjur"} can only be used with models constructed on presence-absence data, the coefficient of determination (\eqn{R^2}{R2}) is calculated by averaging the estimated probability for presences and absences respectively and subtracting the average over the absences from the average over the presences.
#'
#' If \code{type = "efron"} is used, than the \eqn{R^2}{R2} will be calculated using Efron's pseudo-\eqn{R^2}{R2} regardless of the model used. Efron's pseudo-\eqn{R^2}{R2} is calculated as:
#'
#' \deqn{R^2 = 1 - \frac{(y_i-\hat{\pi})^2}{(y_i-\bar{y})^2}{R2=1-((y-pihat)^2)/((y-ybar)^2)}}.
#'
#' where \eqn{\hat{\pi}}{pihat} is the model predicted values calculated on the scale of the response variable.
#'
#' The ordinary least squares ("ols") \eqn{R^2}{R2} and Efron's pseudo-\eqn{R^2}{R2} calculated on Gaussian models lead to the same results.
#'
#' If \code{type = "nakagawa"} is used, than the \eqn{R^2}{R2} will be calculated using Nakagawa and and Schielzeth's \eqn{R^2}{R2} regardless of the model used. Nakagawa and and Schielzeth's \eqn{R^2}{R2} is calculated as:
#'
#' \deqn{R^2 =  \frac{\sigma^2_X + \sum \sigma^2_{Random}}{\sigma^2_X + \sum \sigma^2_{Random} + \sigma^2_{Add} + \sigma^2_{Dist}}{R2 = (varX+sum(varRandom))/(varX+sum(varRandom)+varAdd+varDist)}}.
#'
#' where \eqn{\sigma^2_X}{varX} is the variance associated to the fixed effect (\code{X}), \eqn{\sum\sigma^2_{Random}}{sum(varRandom)} is the sum of the variances associated to the each random effect (\code{Random}), \eqn{\sigma^2_{Add}}{varAdd} is the additive variance  (dispersion) not explained by either the fixed or the random effect and \eqn{\sigma^2_{Dist}}{varDist}. As such, the \eqn{R^2}{R2} implemented in this function is the conditional \eqn{R^2}{R2} of Nakagawa and and Schielzeth's.
#'
#' The calculation of \eqn{\sigma^2_X}{varX} is obtained as follows
#'
#' \deqn{\sigma^2_X = var(X\beta)}{varX = var(X x paramX)}
#'
#' while the calculation of \eqn{\sigma^2_{Random}}{varRandom} is obtained using the same calculation as above but using with latent variables and their associated regression parameters.
#'
#' \deqn{\sigma^2_{Random} = var(H\Lambda)}{varX = var(Latent x paramLatent)}
#'
#' As for \eqn{\sigma^2_{Add}}{varAdd}, it is calculated as follows
#'
#' \deqn{\sigma^2_{Add} = \sum(g(g^{-1}(Y - \hat{Y}))^2)}{varAdd = sum( linkfun(linkinv(Y-Yhat))^2)}
#'
#' where \eqn{g}{linkfun} and \eqn{g^{-1}}{linkinv} are the link and an the inverse link function, respectively.
#'
#' Lastly, \eqn{\sigma^2_{Dist}}{varDist} is link function specific. If the identity link function is used it is 0, if the logit function is used it is \eqn{\pi^2/3}{pi2/3}, if the probit link function is used it is 1 and if the log link function is used it is \eqn{\ln(1/\exp(\beta_0)+1)}{ln(1/exp(beta0)+1)}.
#'
#' The "ols" and "nakagawa" \eqn{R^2}{R2} can also be used to calculate the site's contribution to the \eqn{R^2}{R2}. The argument \code{indepentSite} can be used to decide wether the site's contribution to the \eqn{R^2}{R2} can be calculated. The other types of \eqn{R^2}{R2} cannot be used to calculate the site's contribution to the \eqn{R^2}{R2} because it is either mathematically not possible (Tjur's \eqn{R^2}{R2}) or it gives non-sensical results (Efron's pseudo-\eqn{R^2}{R2} sometimes leads to sites contribution being negative).
#'
#' The adjustement to the coefficient of determination used here is the one proposed by Gelman and Pardoe (2006), which leads to a lower adjusted coefficient of determination than the classical adjustement because it accounts for uncertainty in the model variance. The adjustement used in this function is calculated as follows:
#'
#' \deqn{1 - \frac{n-3}{n-p-2}(1-R^2)}{1-((n-3)/n-p-2)(1-R2)}.
#'
#' @return
#'
#' If \code{averageSp} is \code{TRUE} and the \code{indepentSite} is \code{FALSE}, a single value is returned leading to a community-level coefficient of determination.
#'
#' If \code{averageSp} is \code{FALSE} and the \code{indepentSite} is \code{FALSE} ,a vector with as many values as the number of species is returned leading to specie-level coefficient of determination.
#'
#' If \code{averageSp} is \code{TRUE} and the \code{indepentSite} is \code{TRUE} a vector with as many values as the number of sites is returned presenting the sites contribution to the community-level coefficient of determination.
#'
#' If \code{averageSp} is \code{FALSE} and the \code{indepentSite} is \code{TRUE} a matrix with the same dimension as the community matrix is returned presenting the sites contribution to each species coefficient of determination.
#'
#' Note that when Tjur's or Efron's \eqn{R^2}{R2} are used, \code{indepentSite} will always be considered \code{FALSE}.
#'
#' @references
#'
#' Gelman, A. and I. Pardoe (2006) Bayesian Measures of Explained Variance and Pooling in Multilevel (Hierarchical) Models. \emph{Technometrics} \strong{48}, 241-251.
#'
#' Efron, B. (1978) Regression and ANOVA with Zero-One Data: Measures of Residual Variation., \emph{Journal of the American Statistical Association} \strong{73}, 113-121.
#'
#' Nakagawa, S. and H. Schielzeth (2013) A general and simple method for obtaining \eqn{R^2}{R2} from generalized linear mixed-effects models. \emph{Methods in Ecology and Evolution} \strong{63}, 133–142.
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
Rsquared <- function(hmsc, newdata = NULL, type = c("efron", "ols", "nakagawa", "tjur"), adjust = FALSE, averageSp = TRUE, indepentSite = FALSE) {
  ### Match arguments
  type <- match.arg(type)

  ### Account for newdata
  if (is.null(newdata)) {
    Y <- hmsc$data$Y
  }

	if (!is.null(newdata)) {
	  Y <- newdata$Y
    if(indepentSite){
      stop("The contribution of each site to R2 cannot be calculated with 'newdata'")
    }
  }

  Xvariance <- apply(hmsc$data$X, 2, var)
  if(!any(Xvariance == 0)){
    warning("A model without intercept may lead to a wrong R2")
  }

  ### basic objects
  nsite <- nrow(Y)
  nsp <- ncol(Y)

  ### Calculate model estimates
  if(!(type == "nakagawa")){
    Ypred <- predict(hmsc, newdata=newdata)
  }else{
    Ypred <- predict(hmsc, newdata=newdata, type = "link")
  }
  ### Tjur's R2
  if(type == "tjur"){
    ### A few basic checks
    if(!any(unique(as.vector(Y))%in%c(0,1))){
      stop("Tjur R2 can only be used for binary (0, 1) data")
    }

    if(indepentSite){
      stop("The contribution of each site to R2 cannot be calculated with 'tjur'")
    }

    ### Probit model
    if (any(class(hmsc) == "probit") | any(class(hmsc) == "logit")) {
      Y0 <- Y == 0
      Y1 <- Y == 1

      R2 <- numeric()

      for (i in 1:nsp) {
          R2[i] <- mean(Ypred[Y1[, i], i]) - mean(Ypred[Y0[, i], i])
      }
    }else{
      stop("Tjur's R2 can only calculated on probit and logit models")
    }

    ### Community-level R2
    if (averageSp) {
        R2 <- mean(R2)
    }
  }

  if((type == "ols" & any(class(hmsc) == "gaussian")) | type == "efron"){
    ### Total sums of squares per species
    ssY <- colSums(sweep(Y,2,colMeans(Y),FUN="-")^2)

    ### Residual sums of squares per species
    ssRes <- (Y-Ypred)^2

    ### Result object
    R2 <- matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
    colnames(R2) <- colnames(Y)
    rownames(R2) <- rownames(Y)
    ### Calculate R2
    for (i in 1:nsp) {
      R2[,i] <- 1/nsite - ssRes[,i]/ssY[i]
    }

    ### Global R2
    if(!indepentSite){
      R2 <- colSums(R2)

      ### Community-level R2
      if (averageSp) {
          R2 <- mean(R2)
      }
    }else{
      ### Community-level R2 per site
      if (averageSp) {
        R2 <- rowMeans(R2)
      }
      if(type == "efron"){
        R2 <- colSums(R2)
      }
    }
  }

  if(type == "nakagawa"){
    if(any(class(hmsc) == "gaussian")){
      ### Link and inverse link function
      link <- gaussian(link = "identity")

      ### Distribution specific-variance
      varDist <- 0
    }

    if(any(class(hmsc) == "probit")){
      ### Link and inverse link function
      link <- binomial(link = "probit")

      ### Distribution specific-variance
      varDist <- 1
    }

    if(any(class(hmsc) == "logit")){
      ### Link and inverse link function
      link <- binomial(link = "logit")

      ### Distribution specific-variance
      varDist <- pi^2/3
    }

    if(any(class(hmsc) == "poisson") | any(class(hmsc) == "overPoisson")){
      ### Link and inverse link function
      link <- poisson(link = "log")

      #---------------------------------
      ### Distribution specific-variance
      #---------------------------------
      intercept <- which(Xvariance == 0)[1]

      if(length(intercept) == 0){
        varDist <- rep(log(2), nsp)
      }else{
        paramX0 <- coef(hmsc)$paramX[,intercept]
        varDist <- log(1/exp(paramX0) + 1)
      }
    }

    ### Model variance
    YMeans <- matrix(colMeans(Ypred),
                     nrow = nsite, ncol = nsp, byrow = TRUE)

    varModelSite <- (Ypred - YMeans)^2/(nsite-1)
    varModel <- colSums(varModelSite)

    ### Additive distributional variance
    varAdd <- diag(var(Y-link$linkinv(Ypred)))

    ### Calculate R2
    R2 <- varModelSite / (varModel + varAdd + varDist)

    ### Global R2
    if(!indepentSite){
      R2 <- colSums(varModelSite) / (varModel + varAdd + varDist)

      ### Community-level R2
      if (averageSp) {
          R2 <- mean(R2)
      }
    }else{
      ### Community-level R2 per site
      if (averageSp) {
        R2 <- rowMeans(R2)
      }
    }
  }

  ### Adjustement
  if(adjust){
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

  return(R2)
}
