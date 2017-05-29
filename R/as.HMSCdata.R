#' @title Format data, parameters, and priors to "HMSC" classes
#'
#' @description These functions are used to format the data, parameters, and priors to object of class "HMSC".
#'
#' @param data An object of the class \code{HMSCdata}.
#' @param priors An object of the class \code{HMSCprior}.
#' @param Y A matrix or a data.frame where each column is a species.
#' @param X A matrix or a data.frame where each column is a descriptor of the sites.
#' @param Tr A matrix or a data.frame where each row is a trait characterizing all species.
#' @param Phylo A square correlation or covariance matrix describing the phylogenetic correlation (or covariance) between pairs of all species (see details).
#' @param Auto A data.frame or a list of matrices presenting the spatial (or temporal) coordinates of the samples (see details).
#' @param Random A factor defining a random effect on the samples or a data.frame characterizing multiple random effects.
#' @param scaleX Logical. Whether the columns of X should be centred and divided by the standard deviation. Default is TRUE.
#' @param scaleTr Logical. Whether the rows of Tr should be centred and divided by the standard deviation. Default is TRUE.
#' @param interceptX Logical. Whether a column of 1s should be added to X. Default is TRUE.
#' @param interceptTr Logical. Whether a row of 1s should be added to Tr. Default is TRUE.
#' @param paramX A matrix or data.frame of model parameters defining how species (rows) are characterized by the descriptors (columns).
#' @param meansParamX A vector of parameters defining how an average species reacts to each descriptor in X. In as.HMSCprior, this argument defines the prior information associated to parameters of the same name.
#' @param varMeansParamX In as.HMSCprior, it is a square matrix of parameters defining how \code{meansParamX} varies. This argument is only used when \code{Phylo} is included in the analysis.
#' @param varX Symmetric covariance matrix. Each dimension of this matrix should be equal to the number of explanatory variables.
#' @param paramTr In function as.HMSCparam, it is a matrix of model parameters defining how descriptors (rows) characterizes traits (columns). In function as.HMSCprior, it is a matrix of prior defining how descriptors (rows) characterizes traits (columns).
#' @param varTr In function as.HMSCdata, it is a symmetric covariance matrix. Each dimension of this matrix should be equal to the number of traits. In function as.HMSCprior, it is a symmetric covariance matrix of prior. Each dimension of this matrix should be equal to the number of traits.
#' @param paramPhylo Numeric. In function as.HMSCparam, it defines the importance of phylogeny in structuring the community. In as.HMSCprior, this argument defines prior information associated to the parameter of the same name.
#' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model.
#' @param varDist Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model.
#' @param latent List. Each level of the list includes a matrix of latent (unsampled) variables. There are as many level in the list as there are columns in \code{Random}. The number of rows in each matrix corresponds to the number of levels in each factor in \code{Random}.
#' @param paramLatent List. Each level of the list includes a matrix of parameters for the latent (unsampled) variables. There are as many level in the list as there are columns in \code{Random}. Each matrix has as many rows as there are species in \code{Y} and as many columns as there are latent variables in \code{latent}.
#' @param shrinkLocal Numeric. Hyperparameter of a gamma distribution defining the local shrinkage for latent variables.
#' @param paramShrinkGlobal List. Each level of the list includes a vector of values describing how much shrinkage should be imposed on each latent variables (\code{latent}) as a whole. There are as many level in the list as there are columns in \code{Random}. The number of elements in each vector equals the number of variables in \code{latent}.
#' @param paramAuto Vector. In function as.HMSCparam, it is a vector of numerical values defining the autocorrelated patterns for each factor in \code{Auto}. These values can range from 0 to the largest distance between samples in the autocorrelated level. In as.HMSCprior, this argument defines prior information associated to the parameter of the same name.
#' @param latentAuto List. Each level of the list includes a matrix of autocorrelated latent (unsampled) variables. There are as many level in the list as there are columns in \code{Auto}. The number of rows in each matrix corresponds to the number of levels in each autocorrelated factor in \code{Auto}.
#' @param paramLatentAuto List. Each level of the list includes a matrix of parameters for the autocorrelated latent (unsampled) variables. There are as many level in the list as there are columns in \code{Auto}. Each matrix has as many rows as there are species in \code{Y} and as many columns as there are autocorrelated latent variables in \code{latentAuto}.
#' @param shrinkLocalAuto List. Each level of the list includes a matrix of values describing how much shrinkage should be imposed on each parameter of the autocorrelated latent (unsampled) variables (\code{paramLatentAuto}). There are as many level in the list as there are columns in \code{Auto}. The size of each matrix is the same as the size of \code{paramLatentAuto}.
#' @param paramShrinkGlobalAuto List. Each level of the list includes a vector of values describing how much shrinkage should be imposed on each autocorrelated latent variables (\code{latent}) as a whole. There are as many level in the list as there are columns in \code{Auto}. The number of elements in each vector equals the number of variables in \code{latentAuto}.
#' @param shrinkOverall A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the overall shrinkage of the autocorrelated latent variable is handled (see details).
#' @param shrinkSpeed A vector of length 2 defining the shape and scale hyperparameter of a gamma distribution. The size of the first parameters indicate how the speed of the shrinkage of the autocorrelated latent variable is handled (see details).
#' @param paramAutoWeight Prior. Matrix with a single column defining the importance each value in \code{priorParamAutoDist} may have. The size of \code{priorParamAutoWeight} needs to be equal to the size of \code{priorParamAutoGrid}.
#' @param paramAutoDist Prior. Matrix with a single column defining the potential values given in \code{paramAuto} can take. The size of \code{priorParamAutoGrid} needs to be equal to the size of \code{priorParamAutoWeight}.
#' @param family Character string defining the type of generalized regression model for which priors need to be defined.
#' @param varXDf Numeric. Prior for varX. It is the number of degrees of freedom of the inverse Wishart distribution from which varX is sampled.
#' @param varDistShape Numeric. Prior for the variance of the normal distribution. It is the shape parameter of a gamma distribution from which the variance is sampled.
#' @param varDistScale Numeric. Prior for the variance of the normal distribution. It is the scale parameter of a gamma distribution from which the variance is sampled.
#' @param varXScaleMat Symmetric covariance matrix. Prior for varX. It is the scale matrix of the inverse Wishart distribution from which varX is sampled.
#' @param nsp Numeric. Number of species in \code{Y}.
#'
#' @details
#'
#' Although the structure of the \code{Auto} argument may seem unusual as a list of data.frames, this structures give the possibility to account for different type of autocorrelation structure simultaneously (e.g. spatial and temporal) in the same analysis. Although it is usual for autocorrelation to be accounted for by either 1 or 2 dimensions (so 1 or 2 columns), the analyses carried out in this package can account for autocorrelation at any number of dimensions.
#'
#' Functions \code{as.HMSCparam} and \code{as.HMSCprior} do not need to be used when estimating a model using \code{\link{hmsc}}. In \code{\link{hmsc}}, if \code{as.HMSCparam} and \code{as.HMSCprior} are not defined, the function will define a set of default parameters and priors. \code{as.HMSCparam} and \code{as.HMSCprior} were designed so that some, not necessarily all, priors and parameters need to be defined. For example, if one defines only the prior (parameters) for \code{varTr}, all other priors (parameters) will be defined with default priors (parameter).
#'
#' In \code{as.HMSCparam}, the argument \code{family} is used to define the  \code{paramX} by default using a univariate generalized linear model for each species. Currently, only \code{binomial(link = "probit")} should be used because \code{\link{hmsc}} only estimate probit models.
#'
#' In the hmsc function, the \code{Phylo} argument has to be a square correlation matrix. So, if a covariance matrix is given as \code{Phylo} the \code{as.HMSCdata} function will convert it to a correlation matrix.
#'
#' The choosing the values for \code{shrinkOverall} and \code{shrinkSpeed} is not straightforward. Bhattacharya and Dunson (2011) proposed to use a 2.1 and 1 for \code{shrinkOverall} and 3.1 and 1 for \code{shrinkSpeed}. However, for presence-absence data, it seems that these shrinkage values (especially the scale [first] parameter) are not large enough, which leads to a too many latent variables defined. It is the reason why we suggest to use larger scale values. As such, if they are not predefine, the scale value for \code{shrinkOverall} is 10 whereas the scale value for \code{shrinkSpeed} is 15. These two scale parameters are likely the only one that need to be chosen more carefully to ensure that the model estimation is accurate.
#'
#' \code{shrinkOverall} and \code{shrinkSpeed} are arguments that are used both for the latent variables as well as the autocorrelated latent variables.
#'
#' When predefining the initial value of \code{paramPhylo}, the \code{as.HMSCparam} function may change slightly the value of \code{paramPhylo} to make sure that the \code{hmsc} function runs properly. Since \code{paramPhylo} is estimated using a step function instead of sampling the parameter from a conjugate prior distribution, this step is important to ensure that the code runs well.
#'
#' @return
#'
#' An object of the class \code{HMSCdata} is returned by \code{as.HMSCdata}.
#' An object of the class \code{HMSCparam} is returned by \code{as.HMSCparam}.
#' An object of the class \code{HMSCprior} is returned by \code{as.HMSCprior}.
#'
#' All of these objects present the same output data as the input except for \code{as.HMSCparam} where \code{meansparamX} is also defined.
#'
#' @references Bhattacharya, A. and Dunson, D.B. (2011) Sparse Bayesian infinite factor models. \emph{Biometrika} \strong{98}, 291--306.
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom stats cov2cor
#' @importFrom stats sd
#' @examples
#'
#' #================
#' ### Generate data
#' #================
#' desc <- cbind(1, scale(1:50), scale(1:50)^2)
#' randomEff <- as.factor(rep(letters[1:5], each=10))
#' dataBase <- communitySimul(desc, nsp = 30, Random = randomEff)
#'
#' #=============
#' ### Formatting
#' #=============
#' ### Format data
#' formdata <- as.HMSCdata(Y = dataBase$data$Y, X=desc, Random=randomEff, interceptX = FALSE)
#'
#' ### Format priors
#' formpriors <- as.HMSCprior(formdata, shrinkOverall = c(50, 1), shrinkSpeed = c(20, 1))
#'
#' ### Format parameters
#' formparam <- as.HMSCparam(formdata, formpriors, paramX = dataBase$param$paramX)
#'
#' @keywords manip
#' @keywords classes
#' @export
as.HMSCdata <-
function(Y=NULL, X=NULL, Tr=NULL, Phylo=NULL, Auto=NULL, Random=NULL,
		  formula=NULL,scaleX=TRUE, scaleTr=TRUE, interceptX=TRUE, interceptTr=TRUE){
#### F. Guillaume Blanchet - July 2014, January 2016
##########################################################################################
	Ynames <- colnames(Y)

	if(is.null(Y) & all(is.null(X) & is.null(Random) & is.null(Auto))){
		stop("If 'Y' is NULL, either 'X', 'Random' or 'Auto' can be NULL, not all three of them")
	}

	### Check that neither X or Random are NULL
	if(is.null(X) && is.null(Random) && is.null(Auto)){
		stop("Either 'X', 'Random' or 'Auto' can be NULL, not all three of them")
	}

	### If Auto is a data.frame convert it to a list
	if(!is.null(Auto)){
		if(is.data.frame(Auto)){
			tmp <- Auto
			Auto <- vector("list", length=1)
			Auto[[1]] <- tmp
		}
		if(!is.list(Auto)){
			stop("'Auto' needs to be data.frame or a list")
		}
	}

	### Check for NAs
	if(!is.null(Y)){
		if(any(is.na(Y))){
			stop("There is at least one NA in 'Y'")
		}
	}

	if(!is.null(X)){
		if(any(is.na(X))){
			stop("There is at least one NA in 'X'")
		}
	}

	if(!is.null(Tr)){
		if(any(is.na(Tr))){
			stop("There is at least one NA in 'Tr'")
		}
	}

	if(!is.null(Phylo)){
		if(any(is.na(Phylo))){
			stop("There is at least one NA in 'Phylo'")
		}
	}

	if(!is.null(Auto)){
		if(any(is.na(unlist(Auto)))){
			stop("There is at least one NA in 'Auto'")
		}
	}

	if(!is.null(Random)){
		if(any(is.na(Random))){
			stop("There is at least one NA in 'Random'")
		}
	}

	### Check for non-numeric values
	if(!is.null(Y)){
		if(any(is.na(suppressWarnings(as.numeric(apply(Y, 1, as.character)))))){
			stop("There is at least one non-numeric value in 'Y'")
		}
	}

	if(!is.null(X)){
		if(any(is.na(suppressWarnings(as.numeric(apply(X, 1, as.character)))))){
			stop("There is at least one non-numeric value in 'X'")
		}
	}

	if(!is.null(Tr)){
		if(any(is.na(suppressWarnings(as.numeric(apply(Tr, 1, as.character)))))){
			stop("There is at least one non-numeric value in 'Tr'")
		}
	}

	if(!is.null(Phylo)){
		if(any(is.na(suppressWarnings(as.numeric(apply(Phylo, 1, as.character)))))){
			stop("There is at least one non-numeric value in 'Phylo'")
		}
	}

	#### Check format
	if(!is.null(Y)){
		if(length(dim(Y))!=2){
			stop("'Y' should be a table")
		}
	}

	if(!is.null(X)){
		if(length(dim(X))!=2){
			stop("'X' should be a table")
		}
	}

	if(!is.null(Tr)){
		if(length(dim(Tr))!=2){
			stop("'Tr' should be a table")
		}
	}

	if(!is.null(Phylo)){
		if(nrow(Phylo)!=ncol(Phylo)){
			stop("'Phylo' should be a square table")
		}

		### Check if symmetry
		if(!isSymmetric(Phylo)){
			stop("'Phylo' need to be a symmetric matrix")
		}

		### Check positive definiteness
		if(any(eigen(Phylo)$value<0)){
			stop("'Phylo' need to be a positive definite matrix")
		}

		### Check if Phylo is a covariance matrix, and convert to correlation if it is
		if(any(diag(Phylo)!=1)){
			Phylo <- cov2cor(Phylo)
			print("'Phylo' was converted to a correlation matrix")
		}
	}

	if(!is.null(Auto)){
		if(!all(sapply(Auto, function(x) is.factor(x[, 1])))){
			stop("The first column of the data.frame in each list should be a factor and the other columns should be coordinates")
		}

		### Find the number of columns with coordinates
		AutoCoord <- lapply(lapply(Auto, function(x) x[, -1]), as.matrix)

		if(!all(unlist(lapply(AutoCoord, function(x) apply(x, 2, is.numeric))))){
			stop("When Auto is a list, the first column of the data.frame in each list should be a factor and the other columns should be coordinates")
		}
	}

	if(!is.null(Random)){
		if(is.factor(Random)){
			Random <- data.frame(random=Random)
			print("'Random' was included into a data.frame")
		}else{
			if(is.data.frame(Random)){
				if(!all(mapply(is.factor, Random))){
					stop("If 'Random' is a data.frame, it should only include factors")
				}
			}else{
				stop("'Random' should be a factor or a data.frame")
			}
		}
	}

	#### Check if dimensions of all tables match (This check is carried out only when Y is NULL)
	if(!is.null(Y)){
		if(!is.null(X)){
			if(nrow(Y)!=nrow(X)){
				stop("'X' and 'Y' should have the same number of rows")
			}
		}

		if(!is.null(Tr)){
			if(ncol(Y)!=ncol(Tr)){
				stop("'Y' and 'Tr' should have the same number of columns")
			}
		}

		if(!is.null(Phylo)){
			if(ncol(Y)!=ncol(Phylo)){
				stop("The number of columns of 'Y' should equal the number of rows and columns of 'Phylo'")
			}
		}

		if(!is.null(Auto)){
			nrowAuto <- sapply(Auto, nrow)
			if(any(nrowAuto!=nrow(Y))){
				stop("At least one data.frame of 'Auto' has a different number of rows than 'Y'")
			}
		}

		if(!is.null(Random)){
			if(nrow(Random)!=nrow(Y)){
				stop("'Random' and 'Y' should have a the same number of rows")
			}
		}
	}

	#### Check column names
	if(!is.null(Y)){
		if(is.null(Ynames)){
			colnames(Y) <- paste("y", 1:ncol(Y), sep="")
			print(paste("column names were added to 'Y'"))
		}
	}

	if(!is.null(X)){
		if(is.null(colnames(X))){
			colnames(X) <- paste("x", 1:ncol(X), sep="")
			print("column names were added to 'X'")
		}
	}

	if(!is.null(Tr)){
		if(is.null(colnames(Tr))){
			colnames(Tr) <- paste("y", 1:ncol(Tr), sep="")
			print("column names were added to 'Tr'")
		}
	}

	if(!is.null(Phylo)){
		if(is.null(colnames(Phylo))){
			colnames(Phylo) <- paste("y", 1:ncol(Phylo), sep="")
			print("column names were added to 'Phylo'")
		}
	}

	if(!is.null(Auto)){
		colNamesAuto <- sapply(Auto, function(x) is.null(colnames(x)))
		if(any(colNamesAuto)){
			colnames(Auto[[i]])[1] <- "autoRandom1"
			for(i in which(colNamesAuto)){
				colnames(Auto[[i]])[-1] <- paste("coord", 1:ncol(Auto[[i]])[-1], sep="")
			}
			print("column names were added to at least one section of 'Auto'")
		}
	}

	if(!is.null(colnames(Random))){
		colnames(Random) <- paste("random", 1:ncol(Random), sep="")
		print("column names were added to 'Random'")
	}

	#### Check row names
	if(!is.null(Y)){
		if(is.null(rownames(Y))){
			rownames(Y) <- paste("site", 1:nrow(Y), sep="")
			print(paste("row names were added to 'Y'"))
		}
	}

	if(!is.null(X)){
		if(is.null(rownames(X))){
			rownames(X) <- paste("site", 1:nrow(X), sep="")
			print("row names were added to 'X'")
		}
	}

	if(!is.null(Tr)){
		if(is.null(rownames(Tr))){
			rownames(Tr) <- paste("t", 1:nrow(Tr), sep="")
			print("row names were added to 'Tr'")
		}
	}

	if(!is.null(Phylo)){
		if(is.null(rownames(Phylo))){
			rownames(Phylo) <- paste("y", 1:ncol(Phylo), sep="")
			print("row names were added to 'Phylo'")
		}
	}

	if(!is.null(Auto)){
		rowNamesAuto <- sapply(Auto, function(x) is.null(rownames(x)))
		if(any(rowNamesAuto)){
			for(i in which(rowNamesAuto)){
				rownames(Auto[[i]]) <- paste("site", 1:nrow(Auto[[i]]), sep="")
			}
			print("row names were added to at least one section of 'Auto'")
		}
	}

	if(!is.null(Random)){
		rownames(Random) <- paste("site", 1:nrow(Random), sep="")
	}

	### Check list name
	if(!is.null(Auto)){
		if(is.null(names(Auto))){
			names(Auto) <- paste("auto", 1:length(Auto), sep="")
			print("names were added to each section of 'Auto'")
		}
	}

	### Add an intercept to X and scale
	if(!is.null(X)){
		if(interceptX){
			if(scaleX){
				### Check if any columns have a null variance
				zeroVar <- which(apply(X, 2, sd)==0)
				if(length(zeroVar)!=0){
					X[, -zeroVar] <- scale(X[, -zeroVar])
					warning(paste(colnames(X)[zeroVar], "are explanatory variable(s) with a variance of 0, for this reason no intercept were added, check to make sure this is OK"))
				}else{
					X <- cbind(1, scale(X))
					colnames(X)[1] <- "Intercept"
				}
			}else{
				X <- cbind(1, X)
				colnames(X)[1] <- "Intercept"
			}
		}else{
			if(scaleX){
				### Check if any columns have a null variance
				zeroVar <- which(apply(X, 2, sd)==0)
				if(length(zeroVar)!=0){
					X[, -zeroVar] <- scale(X[, -zeroVar])
					warning(paste(colnames(X)[zeroVar], "are explanatory variable(s) with a variance of 0, make sure this is OK"))
				}else{
					X <- scale(X)
				}
			}
		}
	}

	### Add an intercept to Tr
	if(!is.null(Tr)){
		if(interceptTr){
			if(scaleTr){
				### Check if any columns have a null variance
				zeroVar <- which(apply(Tr, 1, sd)==0)
				if(length(zeroVar)!=0){
					Tr[-zeroVar, ] <- scale(Tr[-zeroVar, ])
					warning(paste(rownames(Tr)[zeroVar], "are trait(s) with a variance of 0, for this reason no intercept were added, check to make sure this is OK"))
				}else{
					Tr <- rbind(1, t(scale(t(Tr))))
					rownames(Tr)[1] <- "Intercept"
				}
			}else{
				Tr <- rbind(1, Tr)
				rownames(Tr)[1] <- "Intercept"
			}
		}else{
			if(scaleTr){
				### Check if any columns have a null variance
				zeroVar <- which(apply(Tr, 1, sd)==0)
				if(length(zeroVar)!=0){
					Tr[-zeroVar, ] <- scale(Tr[-zeroVar, ])
					warning(paste(rownames(Tr)[zeroVar], "are traits with a variance of 0, make sure this is OK"))
				}
			}
		}
	}

	#### Check classes
	if(!is.null(Y)){
		if(!is.matrix(Y)){
			Y <- as.matrix(Y)
			print("'Y' was converted to a matrix")
		}
	}

	if(!is.null(X)){
		if(!is.matrix(X)){
			X <- as.matrix(X)
			print("'X' was converted to a matrix")
		}
	}

	if(!is.null(Tr)){
		if(!is.matrix(Tr)){
			Tr <- as.matrix(Tr)
			print("'Tr' was converted to a matrix")
		}
	}

	if(!is.null(Phylo)){
		if(!is.matrix(Phylo)){
			Phylo <- as.matrix(Phylo)
			print("'Phylo' was converted to a matrix")
		}
	}
	###  ICI ###
	###  ICI ###
	###  ICI ###
	###  ICI ###
	###  ICI ###
	if(!is.null(formula)){
		if()

	}

	#### Return results

	### Including Y
	if(!is.null(Y) & !is.null(X) & is.null(Tr) & is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X)
		attributes(res) <- list(names=c("Y", "X"))
	}

	if(!is.null(Y) & is.null(X) & is.null(Tr) & is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, Random=Random)
		attributes(res) <- list(names=c("Y", "Random"))
	}

	if(!is.null(Y) & is.null(X) & is.null(Tr) & is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, Auto=Auto)
		attributes(res) <- list(names=c("Y", "Auto"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr)
		attributes(res) <- list(names=c("Y", "X", "Tr"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Phylo=Phylo)
		attributes(res) <- list(names=c("Y", "X", "Phylo"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Phylo) & is.null(Tr) & is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Phylo) & is.null(Tr) & !is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Auto=Auto)
		attributes(res) <- list(names=c("Y", "X", "Auto"))
	}

	if(!is.null(Y) & is.null(X) & is.null(Phylo) & is.null(Tr) & !is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("Y", "Auto", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Phylo=Phylo)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Phylo"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Auto=Auto)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Auto"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Phylo=Phylo, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Phylo", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Phylo=Phylo, Auto=Auto)
		attributes(res) <- list(names=c("Y", "X", "Phylo", "Auto"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Tr) & is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Auto", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Phylo=Phylo, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Phylo", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Phylo=Phylo, Auto=Auto)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Phylo", "Auto"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Auto", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Phylo=Phylo, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Phylo", "Auto", "Random"))
	}

	if(!is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(Y=Y, X=X, Tr=Tr, Phylo=Phylo, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("Y", "X", "Tr", "Phylo", "Auto", "Random"))
	}

	### Excluding Y
	if(is.null(Y) & !is.null(X) & is.null(Tr) & is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(X=X)
		attributes(res) <- list(names=c("X"))
	}

	if(is.null(Y) & is.null(X) & is.null(Tr) & is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(Random=Random)
		attributes(res) <- list(names=c("Random"))
	}

	if(is.null(Y) & is.null(X) & is.null(Tr) & is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(Auto=Auto)
		attributes(res) <- list(names=c("Auto"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(X=X, Tr=Tr)
		attributes(res) <- list(names=c("X", "Tr"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(X=X, Phylo=Phylo)
		attributes(res) <- list(names=c("X", "Phylo"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Phylo) & is.null(Tr) & is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Random=Random)
		attributes(res) <- list(names=c("X", "Random"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Phylo) & is.null(Tr) & !is.null(Auto) & is.null(Random)){
		res <- list(X=X, Auto=Auto)
		attributes(res) <- list(names=c("X", "Auto"))
	}

	if(is.null(Y) & is.null(X) & is.null(Phylo) & is.null(Tr) & !is.null(Auto) & !is.null(Random)){
		res <- list(Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("Auto", "Random"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & is.null(Auto) & is.null(Random)){
		res <- list(X=X, Tr=Tr, Phylo=Phylo)
		attributes(res) <- list(names=c("X", "Tr", "Phylo"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Tr=Tr, Random=Random)
		attributes(res) <- list(names=c("X", "Tr", "Random"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(X=X, Tr=Tr, Auto=Auto)
		attributes(res) <- list(names=c("X", "Tr", "Auto"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Phylo=Phylo, Random=Random)
		attributes(res) <- list(names=c("X", "Phylo", "Random"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(X=X, Phylo=Phylo, Auto=Auto)
		attributes(res) <- list(names=c("X", "Phylo", "Auto"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("X", "Auto", "Random"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Tr=Tr, Phylo=Phylo, Random=Random)
		attributes(res) <- list(names=c("X", "Tr", "Phylo", "Random"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & is.null(Random)){
		res <- list(X=X, Tr=Tr, Phylo=Phylo, Auto=Auto)
		attributes(res) <- list(names=c("X", "Tr", "Phylo", "Auto"))
	}

	if(is.null(Y) & !is.null(X) & is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Phylo=Phylo, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("X", "Phylo", "Auto", "Random"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Tr=Tr, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("X", "Tr", "Auto", "Random"))
	}

	if(is.null(Y) & !is.null(X) & !is.null(Tr) & !is.null(Phylo) & !is.null(Auto) & !is.null(Random)){
		res <- list(X=X, Tr=Tr, Phylo=Phylo, Auto=Auto, Random=Random)
		attributes(res) <- list(names=c("X", "Tr", "Phylo", "Auto", "Random"))
	}

	class(res) <- "HMSCdata"
	return(res)
}
