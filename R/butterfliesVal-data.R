#' Butterflies
#'
#' Butterflies of Great Britain
#'
#' @docType data
#'
#' @usage 
#' data(butterfliesTrain)
#' data(butterfliesVal)
#'
#' @details This data was first introduced by Ovaskainen et al. (2016) and reanalyzed by Ovaskainen et al. (submitted). In all, 55 butterflies species were sampled on a 10 km by 10 km resolution grid across Great Britain. In all 2609 grid cells were considered, which were divided into 300 training data (\code{butterfliesTrain}) and 2309 validation data (\code{butterfliesVal}).
#'
#' @format Both the training data (\code{butterfliesTrain}) and the validation data (\code{butterfliesVal}) are objects of class \code{HMSCdata} that include
#' \itemize{
#'     \item{\code{Y}}{ A community \code{matrix} including 55 butterfly species with 300 training data (\code{butterfliesTrain}) and 2309 validation data (\code{butterfliesVal})}
#'     \item{\code{X}}{ A \code{matrix} of 5 covariates (intercept + (1) the number of growing degree days above 5 Celsius and the percentage of the grid cell cover that consists of (2) broadleaved woodland, (3) coniferous woodland and (4) calcareous substrates)}
#'     \item{\code{Tr}}{ A \code{matrix} of 3 traits (Whether the butterfly species is considered to be found across Great Britain, is a specieslist species or is a migrant species)}
#'     \item{\code{Phylo}}{ A \code{matrix} describing the phylogenetic correlations between all 55 species}
#'     \item{\code{Auto}}{A \code{data.frame} that includes the x and y spatial coordinates of each grid cell}
#' }
#'
#' @references
#' Ovaskainen, O., D. B. Roy, R. Fox and B. J. Anderson. 2016. Uncovering hidden spatial structure in species communities with spatially explicit joint species distribution models. \emph{Methods in Ecology and Evolution} \strong{7}:428-436.
#' 
#' Ovaskainen, O., G. Tikhonov, A. Norberg, F. G. Blanchet, L. Duan, D. Dunson, T. Roslin and N. Abrego. Submitted. How to make more out of community ecology data? A road map and its implementation as models and software
#'
#' @keywords datasets
"butterfliesVal"
