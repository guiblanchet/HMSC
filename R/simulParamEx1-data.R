#' Simulated data and parameters
#'
#' Simulated data and the parameters used to construct the simulated data.  These data and their associated parameters were constructed to illustrate some of the functionality of the R package.
#'
#' @docType data
#'
#' @usage 
#' data(simulEx2)
#' data(simulParamEx2)
#'
#' @format 
#' \code{simulEx2} is an object of class \code{HMSCdata} that include:
#' \itemize{
#'     \item{\code{Y}}{ A community \code{matrix}}
#'     \item{\code{X}}{ A \code{matrix} of three explanatory (environmental) variables}
#'     \item{\code{Random}}{ A \code{data.frame} that includes two factors presenting a site-level random effect and a 5-level random effect}
#' }
#' \code{simulParamEx2} is an object of class \code{HMSCparam} that include:
#' \itemize{
#'     \item{\code{paramX}}{ A \code{matrix} that contains the parameters of each covariates (\code{X} in \code{simulEx1}) associated to each simulated species. }
#'     \item{\code{meansParamX}}{ A \code{matrix} that contains a parameter associated to each covariate(\code{X} in \code{simulEx1}). }
#'     \item{\code{varX}}{ A covariance \code{matrix}. }
#'     \item{\code{precX}}{ A precision \code{matrix}, which is the inverse of \code{varX}. Technically this matrix is redundant with \code{varX} but because it is used often in the parameter estimation, it is convenient for \code{\link{as.HMSCparam}} to construct it. }
#'     \item{\code{latent}}{ A \code{list} that contain 2 matrices of two latent variables. }
#'     \item{\code{paramLatent}}{ A \code{list} that contain 2 matrices, which include the parameters associated to the latent variables in \code{latent}. }
#' }
#'
#' @keywords datasets
"simulParamEx1"