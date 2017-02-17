#' Bryophytes
#'
#' Epiphytic bryophytes growing on aspen trees in central Finland
#'
#' @docType data
#'
#' @usage data(bryophytes)
#'
#' @details This data was first introduced by Olden et al. (2014) and reanalyzed by Ovaskainen et al. (submitted). In all, 60 bryophyte species were sampled at 204 sites in 14 natural and 14 harvested forests in central Finland. 
#'
#' @format An object of class \code{HMSCdata} that include
#' \itemize{
#'     \item{\code{Y}}{ A community \code{matrix} including 60 bryophytes species sampled at 204 sites}
#'     \item{\code{X}}{ A \code{matrix} of 5 covariates (intercept + 4 covariates describing the size of the logs on which bryophytes were sampled and the type of forests)}
#'     \item{\code{Tr}}{ A \code{matrix} of 6 traits (intercept + 5 life forms the different species can take)}
#'     \item{\code{Phylo}}{ A square symmetric correlation matrix describing the phylogenetic correlations between each pairs of species.}
#'     \item{\code{Random}}{ A \code{data.frame} that includes two factors presenting a sample-level random effect, a 28 levels random effect describing the forest for each site.}
#' }
#'
#' @note In this data, the covariates (\code{X}) have not been scaled. Because some of the analyses in the manual require that the covariates be scaled, use the following code to perform the analysis.
#'
#' \code{data(bryophytes)}
#' 
#' \code{bryophyteNoScaleX <- bryophytes}
#'
#' \code{bryophytes <- as.HMSCdata(Y=bryophytes$Y, X = bryophytes$X, Tr = bryophytes$Tr,
#'                                 Phylo = bryophytes$Phylo, Random = bryophytes$Random,
#'                                 scaleX = TRUE, scaleTr = FALSE,
#'                                 interceptX = FALSE, interceptTr = FALSE)}
#'
#' @references
#' Olden, A., O. Ovaskainen, J. S. Kotiaho, S. Laaka-Lindberg and P. Halme. 2014. Bryophyte species richness on retention aspens recovers in time but community structure does not. \emph{PLoS ONE} \strong{9}:e93786.
#' 
#' Ovaskainen, O., G. Tikhonov, A. Norberg, F. G. Blanchet, L. Duan, D. Dunson, T. Roslin and N. Abrego. Submitted. How to make more out of community ecology data? A road map and its implementation as models and software
#'
#' @keywords datasets
"bryophytes"