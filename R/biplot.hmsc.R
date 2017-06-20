#' @title Biplot for HMSC
#'
#' @description Draw a biplot for HMSC analysis carried out on a model with that includes at least one random effect
#'
#' @param x An object of the class \code{hmsc}.
#' @param Random Random effect for which the biplot needs to be drawn. This could be either a character string defining the name of the random effect to consider or a value defining the order of the random effect to consider.
#' @param choice Latent variables to show.
#' @param display Scores shown. These must some of the alternatives "species" for species scores, and/or "sites" for site scores.
#' @param type Type of plot: partial match to text for text labels, points for points, and none for setting frames only. If omitted, text is selected for smaller data sets, and points for larger. Can be of length 2 (e.g. type = c("text", "points")), in which case the first element describes how species scores are handled, and the second how site scores are drawn.
#' @param xlim The x limits (min, max) of the plot.
#' @param ylim The y limits (min, max) of the plot.
#' @param col Colours used for sites and species (in this order). If only one colour is given, it is used for both.
#' @param \dots Other parameters for plotting functions.
#'
#' @details
#'
#' In the current biplot, the "species" scores are drawn as points and the point position is obtained from \code{paramLatent}. Similarly, the "sites" scores are also point and their location is obtained from \code{latent}.
#'
#' All the calculations are carried out on the estimation part of the model, the burning part is discarded.
#'
#' Also, the current function is not yet design to draw biplot using any other random effect than the non-autocorrelated random effect.
#'
#' Graphical parameters can also be given to biplot: the size of xlabs and ylabs is controlled by cex.
#'
#' @return
#'
#' An ordination plot for the unconstrained part of the HMSC analysis
#'
#' @references
#'
#' Hui, F. K. C., S. Taskinen, S. Pledger, S. D. Foster, and D. I. Warton. 2015. Model-based approaches to unconstrained ordination. \emph{Methods in Ecology and Evolution} \strong{6}, 399–411.
#'
#' Warton, D. I., F. G. Blanchet, R. B. O’Hara, O. Ovaskainen, S. Taskinen, S. C. Walker, and F. K. C. Hui. 2015. So Many Variables: Joint Modeling in Community Ecology. \emph{Trends in Ecology & Evolution} \strong{30}:766–779.
#'
#' @author F. Guillaume Blanchet
#'
#' @examples
#'
#' #==================
#' ### Simulating data
#' #==================
#' desc <- matrix(1, nrow = 50, ncol = 1)
#' random1 <- as.factor(1:50)
#' random2 <- as.factor(rep(letters[1:2], each = 25))
#' randEff <- data.frame(rand1 = random1, rand2 = random2)
#' nspecies <- 10
#'
#' comm <- communitySimul(X = desc, Random = randEff, nsp = nspecies)
#'
#' #=================================
#' ### Formatting data and parameters
#' #=================================
#' formdata <- as.HMSCdata(Y = comm$data$Y, X = desc, Random = randEff, interceptX = FALSE)
#'
#' #==============
#' ### Build model
#' #==============
#' model <- hmsc(formdata, niter = 200, nburn = 100, thin = 1, verbose = FALSE)
#'
#' #===============================
#' ### Construct correlation matrix
#' #===============================
#' biplot(model,Random=1)
#'
#' @keywords hplot
#' @export
biplot.hmsc <-
function(x, Random, choice = c(1, 2), display = c("sites", "species"), type, xlim, ylim, xlab, ylab, col = c(1, 2) , ...){
#### F. Guillaume Blanchet
##########################################################################################
	### General check
	if(!inherits(x, "hmsc")){
		stop("'biplot.hmsc' is only for objects of class 'hmsc'")
	}

	if(missing(Random)){
		stop("'Random' need to be either a character string or a numerical value")
	}

	if(length(Random)!=1){
		stop("'Random' needs to have a length of 1")
	}

	if(is.character(Random)){
		Random <- which(colnames(x$data$Random)==Random)
		if(length(Random)==0){
			stop("'Random' is not a column name of 'x$data$Random'")
		}
	}

	if(length(choice)!=2){
		stop("'choice' shoulf have a length of 2")
	}

	if(!any(names(x$results$estimation)=="paramLatent" | names(x$results$estimation)=="paramLatentAuto")){
		stop("Only works if 'paramLatent' was estimated (if there were random effects included in the  HMSC model)")
	}

	### Organize arguments
	TYPES <- c("text", "points", "none")

	display <- match.arg(display, several.ok = TRUE)

	if (length(col) == 1){
		col <- c(col,col)
	}

	### Basic objects
	nlevelRandom<-nlevels(x$data$Random[,Random])
	nsp<-ncol(x$data$Y)

	### Number of iterations performed for the model estimation
	nest<-nrow(x$results$estimation$paramLatent)

	### Find the number of latent variables constructed for all iterations
	nLatent<-numeric()
	for(i in 1:nest){
		nLatent[i]<-ncol(x$results$estimation$paramLatent[[i,Random]])
	}
	maxLatent<-max(nLatent)

	if(max(choice)>maxLatent){
		stop("The largest latent variable to show (in choice) is too large ")
	}

	### Organize the latent variables and their associated parameters
	latent<-array(0,dim=c(nlevelRandom,2,nest))
	paramLatent<-array(0,dim=c(nsp,2,nest))

	### This require that the matrices be rotated in case there are some sign issues
	for(i in 1:nest){
		rotationMat<-orthProcrustesRotMat(x$results$estimation$latent[[i,Random]][,choice],x$results$estimation$latent[[nest,Random]][,choice])

		latent[,,i]<-x$results$estimation$latent[[i,Random]][,choice]%*%rotationMat
		paramLatent[,,i]<-x$results$estimation$paramLatent[[i,Random]][,choice]%*%rotationMat
	}

	### Average over latent and paramLatent
	latentMean<-apply(latent,1:2,mean)
	colnames(latentMean)<-paste("latent",choice)
	rownames(latentMean)<-levels(x$data$Random[,Random])

	paramLatentMean<-apply(paramLatent,1:2,mean)
	colnames(paramLatentMean)<-paste("latent",choice)
	rownames(paramLatentMean)<-colnames(x$data$Y)

	### If type is not given
	if (missing(type)) {
      nitlimit <- 80
      nit <- max(nrow(latentMean), nrow(paramLatentMean))
      if (nit > nitlimit){
				type <- rep("points", 2)
			}else{
				type <- rep("text", 2)
			}
  }else{
		type <- match.arg(type, TYPES, several.ok = TRUE)
	}
	if(length(type) < 2){
  	type <- rep(type, 2)
	}

	### Define the range of the plot
	if(missing(xlim)){
    xlim <- range(latentMean[, 1], paramLatentMean[, 1], na.rm = TRUE)
	}
	if(missing(ylim)){
    ylim <- range(latentMean[, 2], paramLatentMean[, 2], na.rm = TRUE)
	}

	### Define the labels
	if(missing(xlab)){
    xlab <- paste("LV",choice[1])
	}
	if(missing(ylab)){
		ylab <- paste("LV",choice[2])
	}

	### Draw the basis of the plot
	plot(0,0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "n", asp = 1, ...)
	abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)

	### Draw sites
	if (any(display == "sites")){
		if(type[1] == "points"){
			points(latentMean,col=col[1])
		}
    if(type[1] == "text"){
			text(latentMean, rownames(latentMean),col = col[1], cex = 0.7)
		}
  }

	### Draw species
	if(any(display == "species")){
    if(type[2] == "points"){
			points(paramLatentMean,col=col[2])
		}
    if(type[2] == "text"){
			text(paramLatentMean, rownames(paramLatentMean),col = col[2], cex = 0.7)
		}
  }

	res<-list(sites=latentMean,species=paramLatentMean)

	class(res) <- "ordiplot"
  invisible(res)
}
