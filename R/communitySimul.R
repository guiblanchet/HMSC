#' @title Simulate a species community
#'
#' @description Simulates a species community assuming a linear model with the proposed descriptors, the species traits and one or more random effects.
#'
#' @param X Matrix of descriptors defining the species. Rows are sites and the columns are descriptors.
#' @param Tr Matrix defining species traits. Rows are traits and columns are species.
#' @param Phylo Square correlation matrix defining the phylogenetic relationship among pairs of species.
#' @param Random A factor or a data.frame that includes many factors defining a random effect on the sites.
#' @param Auto A data.frame that includes a factor as its first column and spatial/temporal coordinates for the other columns. If multiple autocorrelated random effects are considered a list is given where each element of the list is a data.frame as defined previously.
#' @param nsp Numeric. Number of species to be simulated in the community.
#' @param family A character string defining the error distribution and link function to be used when simulating the data. Only "gaussian","probit", "poisson" and "nbinomial" have been implemented so far.
#' @param paramDist Numeric. A value defining the standard deviation (if \code{family="gaussian"}) or the dispersal (if  \code{family="nbinomial"}) matrix of model parameters defining how each species (rows) is characterized by each descriptors (columns).
#' @param paramX A matrix of model parameters defining how each species (rows) is characterized by each descriptors (columns).
#' @param paramTr A matrix of model parameters defining how each descriptors (rows) characterizes species traits (columns).
#' @param paramPhylo A numeric value defining the importance of phylogeny in structuring species relationship with the environment.
#' @param meansParamX Vector of means used to generate the model parameters. (See details)
#' @param varX Covariance matrix used to generate the model parameters. (See details)
#' @param latent A list of matrices where each of set defines a random effect on the sites (rows) characterized by different latent variables (columns).
#' @param paramLatent A list of matrices where each of set defines the parameters of a random effect on the sites (rows) characterized by different latent variables (columns).
#' @param spCor A species by species correlation matrix defining the correlations among species in the community.
#' @param shrinkLocal A list of matrices defining local shrinkage parameters for each parameter of the latent variables associated to each random effect.
#' @param paramShrinkGlobal A list of vectors defining the independent global shrinkage parameters for each latent variable.
#' @param latentAuto A list of matrices where each of set defines an autocorrelated random effect on the sites (rows) characterized by different autocorrelated latent variables (columns).
#' @param paramLatentAuto A list of matrices where each of set defines the parameters of an autocorrelated random effect on the sites (rows) characterized by different autocorrelated latent variables (columns).
#' @param paramAuto A list of numerical values defining the importance of the autocorrelation for each autocorrelated latent variables \code{latentAuto}. These values can range from 0 to the largest distance between samples in the autocorrelated level. (See details)
#' @param shrinkLocalAuto A list of matrices defining local shrinkage parameters for each parameter of the autocorrelated latent variables associated to each autocorrelated random effect.
#' @param paramShrinkGlobalAuto A list of vectors defining the independent global shrinkage parameters for each autocorrelated latent variable.
#'
#' @details
#'
#' This function can be used to simulate a community from randomly proposed or fixed parameters.
#'
#' The descriptors in \code{X} are used without any modifications (or additions) to simulate the species community. As such, a column of 1 should be included in \code{X} for the model used to construct the community to include an intercept.
#'
#' The values in \code{meansParamX} and \code{varX} are used as parameter of a multivariate normal distribution to generate the model's parameters (\code{\link[MASS]{mvrnorm}} is used in the function). When \code{paramX} is set to \code{NULL}, the \code{meansParamX} and the \code{varX} will be randomly generated if they are also set to \code{NULL}. When values are given to \code{paramX} the values of the \code{meansParamX} and the \code{varX} are not used and if either is set to \code{NULL}, no data will be generated for either set of parameter. When generated, the values for \code{meansParamX} are randomly sampled from \code{\link{rnorm}} and the values for \code{varX} are randomly sampled from an inverse (\code{solve})  Wishart distribution (\code{\link{rWishart}}).
#'
#' Note that \code{meansParamX} can be calculated directly by \code{Tr} with \code{paramTr}. As such, \code{meansParamX} is made available as an extension. If it is given to the function but \code{paramTr} is not, than \code{paramTr} is calculated from \code{meansParamX} and vice versa. If both \code{meansParamX} and \code{paramTr} are given to the function and there is a mismatch in the parameters calculated, \code{meansParamX} will take precedent.
#'
#' All the parameters associated to the autocorrelated random effect (\code{Auto},\code{latentAuto},\code{paramLantentAuto} and \code{paramAuto}) use the distance between groups (levels of the factors in \code{Auto}) to define the parameters and latent variables. When multiple samples are within a group, the coordinates related to the different samples are averaged both calculating the distances between this group and other groups.
#'
#' Currently, this function simulates four types of community:
#'
#' \itemize{
#' 	\item{\code{gaussian}}{ For normally distributed data.}
#' 	\item{\code{probit}}{ For presence/absence of species using a probit link.}
#' 	\item{\code{poisson}}{ For count data of species using a log link.}
#' 	\item{\code{nbinomial}}{ For count data with many zeros.}
#' }

#' @return
#'
#' The functions \code{communitySimulH} and \code{communitySimulHT} return an object of class \code{communitySimul} with the following components:
#' \item{data}{an object of class HMSCdata}
#' \item{param}{an object of class HMSCparam}
#' \item{sd}{The standard deviation used to simulate normally distributed data. This value is the same as paramDist and appears only when the data is simulated using \code{gaussian}}
#' \item{size}{The dispersal parameter used when simulating data distributed following a negative binomial distribution. This value is the same as paramDist and appears only when the data is simulated using \code{nbinomial}}
#' \item{probMat}{A matrix that defines the occurrence probability of each species for each site.}
#'
#' @author F. Guillaume Blanchet
#'
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats rWishart
#' @importFrom stats rnbinom
#' @importFrom stats rpois
#' @importFrom stats pnorm
#' @examples
#'
#' ### Construct some descriptors
#' #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#' ### Simulate community from random parameters
#' #=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
#' desc <- cbind(1, scale(1:50), scale(1:50)^2)
#' traits <- rbind(1, 1:20/5)
#' random1 <- as.factor(1:50)
#' random2 <- as.factor(rep(letters[1:2], each = 25))
#' randEff <- data.frame(rand1 = random1, rand2 = random2)
#'
#' ### Simulate presence-absence community data
#' commDescTraitRandEffProbit <- communitySimul(X = desc, Tr = traits, Random = randEff, nsp = 20)
#'
#' #--------------------
#' ### Fixing parameters
#' #--------------------
#' ## ProbMat does not change
#' TrueParamX <- commDescTraitRandEffProbit$param$paramX
#' commParamX <- communitySimul(X = desc, paramX = TrueParamX)
#'
#' ### ProbMat changes
#' TrueParamTr <- commDescTraitRandEffProbit$param$paramTr
#' commParamTr <- communitySimul(X = desc, Tr = traits, paramTr = TrueParamTr)
#'
#' ### ProbMat changes
#' commParamLatent <- communitySimul(X = desc, Random = randEff, paramLatent = commDescTraitRandEffProbit$param$paramLatent, latent = commDescTraitRandEffProbit$param$latent)
#'
#' ### ProbMat changes
#' speciesCor <- cov2cor(solve(rWishart(1, 22, diag(20))[, , 1]))
#' commParamSpCor <- communitySimul(X = desc, spCor = speciesCor)
#'
#' @keywords datagen
#' @keywords htest
#' @importFrom MASS mvrnorm
#' @importFrom stats rgamma
#' @importFrom stats dist
#' @export
communitySimul<-function(X=NULL,Tr=NULL,Phylo=NULL,Random=NULL,Auto=NULL,nsp=NULL,family="probit",paramDist=NULL,paramX=NULL,paramTr=NULL,paramPhylo=NULL,meansParamX=NULL,varX=NULL,latent=NULL,paramLatent=NULL,spCor=NULL,shrinkLocal=NULL,paramShrinkGlobal=NULL,latentAuto=NULL,paramLatentAuto=NULL,paramAuto=NULL,shrinkLocalAuto=NULL,paramShrinkGlobalAuto=NULL){
#### F. Guillaume Blanchet - July 2014, October 2014, January 2016, May 2016
##########################################################################################
	### This makes the outlierSp inactive
	outlierSp<-NULL
	outlierWeight<-0.001

	#================================================
	### Standard format of the X matrix if X is given
	#================================================
	if(!is.null(X)){
		X<-as.matrix(X)
		nsite<-nrow(X)
		nparamX<-ncol(X)

		if(is.null(colnames(X))){
			colnames(X)<-paste("x",1:nparamX,sep="")
		}

		if(is.null(rownames(X))){
			rownames(X)<-paste("site",1:nsite,sep="")
		}

		#=====================
		### If paramX is given
		#=====================
		if(!is.null(paramX)){
			if(is.vector(paramX)){
				nsp<-length(paramX)
				paramX<-matrix(paramX,ncol=1)
			}else{
				nsp<-nrow(paramX)
			}
		}
	}

	#======================
	### If paramTr is given
	#======================
	if(!is.null(paramTr)){
		if(is.null(Tr)){
			stop("'Tr' needs to be given if 'paramTr' is given")
		}
		nsp<-ncol(Tr)
	}

	#=========================
	### If paramPhylo is given
	#=========================
	if(!is.null(paramPhylo)){
		if(is.null(Phylo)){
			stop("'Phylo' needs to be given if 'paramPhylo' is given")
		}
		nsp<-ncol(Phylo)
	}

	#==========================
	### If paramLatent is given
	#==========================
	if(!is.null(paramLatent)){
		if(is.null(latent) && is.null(Random)){
			stop("If 'paramlatent' is given, 'latent' and 'Random' also needs to be given")
		}
		nsp<-nrow(paramLatent[[1]])
	}

	#=====================
	### If latent is given
	#=====================
	if(!is.null(latent)){
		if(is.null(paramLatent) && is.null(Random)){
			stop("If 'latent' is given, 'paramLatent' and 'Random' also needs to be given")
		}
	}

	#====================
	### If spCor is given
	#====================
	if(!is.null(spCor)){
		if(!is.null(latent) && !is.null(paramLatent)){
			stop("If 'spCor' is given, 'latent' and 'paramLatent' do not need to be given")
		}
		if(!is.null(Random)){
			stop("If 'spCor' is given, 'Random' does not need to be given")
		}
		nsp<-nrow(spCor)

		### Basic check
		if(min(spCor) < -1 | max(spCor) > 1){
			stop("'spCor' is a correlation matrix, so that values in it should range from -1 to 1")
		}

		if(nrow(spCor)!=ncol(spCor)){
			stop("'spCor' needs to be a correlation matrix with as many rows as columns")
		}

		if(!is.matrix(spCor)){
			spCor<-as.matrix(spCor)
		}

		### Add row and column names to spCor
		if(is.null(colnames(spCor))){
			colnames(spCor)<-paste("sp",1:nsp,sep="")
		}

		if(is.null(rownames(spCor))){
			rownames(spCor)<-paste("sp",1:nsp,sep="")
		}
	}

	#==========================
	### If shrinkLocal is given
	#==========================
	if(!is.null(shrinkLocal)){
		if((is.null(paramLatent) & is.null(latent)) && is.null(Random)){
			stop("If 'shrinkLocal' is given, 'paramLatent', 'latent' and 'Random' also needs to be given")
		}
	}

	#================================
	### If paramShrinkGlobal is given
	#================================
	if(!is.null(paramShrinkGlobal)){
		if((is.null(paramLatent) & is.null(latent)) && is.null(Random)){
			stop("If 'paramShrinkGlobal' is given, 'paramLatent', 'latent' and 'Random' also needs to be given")
		}
	}

	#========================
	### If paramAuto is given
	#========================
	if(!is.null(paramAuto)){
		if(is.null(latentAuto) && is.null(paramLatentAuto) && is.null(Auto)){
			stop("If 'paramAuto' is given 'latentAuto' and 'paramLatentAuto' also needs to be given")
		}

		nparamAuto<-lapply(paramAuto,length)
		nparamLatentAuto<-lapply(paramLatentAuto,ncol)
		nLatentAuto<-lapply(latentAuto,ncol)

		if(!((length(nparamAuto)==length(nparamLatentAuto))==(length(nLatentAuto)==length(nparamAuto)))){
			stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same length")
		}else{
			if(!(nparamAuto%in%nparamLatentAuto)==(nparamAuto%in%nLatentAuto)){
				stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same values")
			}
		}
	}

	#==============================
	### If paramLatentAuto is given
	#==============================
	if(!is.null(paramLatentAuto)){
		if(is.null(latentAuto) && is.null(paramAuto) && is.null(Auto)){
			stop("If 'paramLatentAuto' is given, 'paramAuto', 'latentAuto' and 'Auto' also needs to be given")
		}
		nsp<-nrow(paramLatentAuto[[1]])

		nparamAuto<-lapply(paramAuto,length)
		nparamLatentAuto<-lapply(paramLatentAuto,ncol)
		nLatentAuto<-lapply(latentAuto,ncol)

		if(!((length(nparamAuto)==length(nparamLatentAuto))==(length(nLatentAuto)==length(nparamAuto)))){
			stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same length")
		}else{
			if(!(nparamAuto%in%nparamLatentAuto)==(nparamAuto%in%nLatentAuto)){
				stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same values")
			}
		}
	}

	#=========================
	### If latentAuto is given
	#=========================
	if(!is.null(latentAuto)){
		if(is.null(paramLatentAuto) && is.null(paramAuto) && is.null(Auto)){
			stop("If 'latentAuto' is given, 'paramLatentAuto', 'paramAuto' and 'Auto' also needs to be given")
		}

		nparamAuto<-lapply(paramAuto,length)
		nparamLatentAuto<-lapply(paramLatentAuto,ncol)
		nLatentAuto<-lapply(latentAuto,ncol)

		if(!((length(nparamAuto)==length(nparamLatentAuto))==(length(nLatentAuto)==length(nparamAuto)))){
			stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same length")
		}else{
			if(!(nparamAuto%in%nparamLatentAuto)==(nparamAuto%in%nLatentAuto)){
				stop("'paramAuto','latentAuto' and 'paramLatentAuto' should all have the same values")
			}
		}
	}

	#==============================
	### If shrinkLocalAuto is given
	#==============================
	if(!is.null(shrinkLocalAuto)){
		if((is.null(paramLatentAuto) & is.null(latentAuto)) && is.null(Auto)){
			stop("If 'shrinkLocalAuto' is given, 'paramLatentAuto', 'latentAUto' and 'Auto' also needs to be given")
		}
	}

	#====================================
	### If paramShrinkGlobalAuto is given
	#====================================
	if(!is.null(paramShrinkGlobalAuto)){
		if((is.null(paramLatentAuto) & is.null(latentAuto)) && is.null(Auto)){
			stop("If 'paramShrinkGlobalAuto' is given, 'paramLatentAuto', 'latentAUto' and 'Auto' also needs to be given")
		}
	}

	#============================================
	### If X, latent and latentAuto are both NULL
	#============================================
	if(is.null(X) & is.null(latent) & is.null(Random) & is.null(Auto) & is.null(latentAuto)){
		stop("At least one of 'X', 'Random', 'latent', 'Auto' or 'latentAuto' needs to be given")
	}

	#===================================
	### Standard format of the Tr matrix
	#===================================
	if(!is.null(Tr)){
		if(!is.vector(Tr)){
			Tr<-as.matrix(Tr)
		}else{
			Tr<-matrix(Tr,nrow=1)
		}
		### Define nsp if it is not yet defined
		if(!is.null(nsp)){
			if(ncol(Tr)!=nsp){
				stop("'nsp' should be the same as the number of column of Tr")
			}
		}else{
			nsp<-ncol(Tr)
		}

		### Organize Tr
		nTr<-nrow(Tr)

		if(is.null(colnames(Tr))){
			colnames(Tr)<-paste("sp",1:nsp,sep="")
		}
		if(is.null(rownames(Tr))){
			rownames(Tr)<-paste("t",1:nTr,sep="")
		}
	}

	#======================================
	### Standard format of the Phylo matrix
	#======================================
	if(!is.null(Phylo)){
		if(ncol(Phylo)!=nrow(Phylo)){
			stop("'Phylo' should be a square correlation matrix")
		}

		### Check if symmetry
		if(!isSymmetric(Phylo)){
			stop("'Phylo' need to be a symmetric matrix")
		}

		### Check positive definiteness
		if(any(eigen(Phylo)$value<0)){
			stop("'Phylo' need to be a positive definite matrix")
		}

		### Check if Phylo is a covariance matrix and convert to correlation if it is
		if(any(diag(Phylo)!=1)){
			stop("'Phylo' need to be a correlation matrix")
		}

		### Define nsp if it is not yet defined
		if(!is.null(nsp)){
			if(ncol(Phylo)!=nsp){
				stop("'nsp' should be the same as the number of column (and rows) of 'Phylo'")
			}
		}else{
			nsp<-ncol(Phylo)
		}

		if(is.null(colnames(Phylo))){
			colnames(Phylo)<-paste("sp",1:nsp,sep="")
		}
		if(is.null(rownames(Phylo))){
			rownames(Phylo)<-paste("sp",1:nsp,sep="")
		}
	}

	#=======================================
	### Standard format of the Random effect
	#=======================================
	if(!is.null(Random)){
		if(is.factor(Random)){
			Random<-data.frame(random1=Random)
			rownames(Random)<-paste("site",1:nsite,sep="")
			### Number or random effects
			nRandom<-1
		}else{
			if(is.data.frame(Random)){
				if(!all(mapply(is.factor,Random))){
					stop("If 'Random' is a data.frame, it should only include factors")
				}
				### Number or random effects
				nRandom<-ncol(Random)

				colnames(Random)<-paste("random",1:nRandom,sep="")
				rownames(Random)<-paste("site",1:nsite,sep="")
			}else{
				stop("'Random' should be a factor or a data.frame")
			}
		}
		nsite<-nrow(Random)
	}

	#======================================================
	### Standard format of the autocorrelated Random effect
	#======================================================
	if(!is.null(Auto)){
		if(is.data.frame(Auto)){
			nsite<-nrow(Auto)

			if(!is.factor(Auto[,1])){
				stop("When Auto is a data.frame, its first column should be a factor and the other columns should be coordinates")
			}
			if(!all(apply(as.matrix(Auto[,-1]),2,is.numeric))){
				stop("When Auto is a data.frame, its first column should be a factor and the other columns should be coordinates")
			}

			### Name stuff
			rownames(Auto)<-paste("site",1:nsite,sep="")
			colnames(Auto)[1]<-paste("autoRandom",1)
			colnames(Auto)[-1]<-paste("coord",1:(ncol(Auto)-1))

			Auto<-list(auto1=Auto)
			### Number or autocorrelated random effects
			nAuto<-1
		}else{
			if(is.list(Auto)){
				nsite<-nrow(Auto[[1]])

				if(!all(mapply(is.data.frame,Auto))){
					stop("If 'Auto' is a list, it should only include data.frames")
				}

				if(!all(sapply(Auto,function(x) is.factor(x[,1])))){
					stop("When Auto is a list, the first column of the data.frame in each list should be a factor and the other columns should be coordinates")
				}

				### Find the number of columns with coordinates
				AutoCoord<-lapply(lapply(Auto,function(x) x[,-1]),as.matrix)

				if(!all(unlist(lapply(AutoCoord,function(x) apply(x,2,is.numeric))))){
					stop("When Auto is a list, the first column of the data.frame in each list should be a factor and the other columns should be coordinates")
				}

				### Number or autocorrelated random effects
				nAuto<-length(Auto)

				### Name stuff
				for(i in 1:nAuto){
					rownames(Auto[[i]])<-paste("site",1:nsite,sep="")
					colnames(Auto[[i]])[1]<-"autoRandom1"
					colnames(Auto[[i]])[-1]<-paste("coord",1:(ncol(Auto[[i]])-1),sep="")
				}

				names(Auto)<-paste("auto",1:nAuto,sep="")

			}else{
				stop("'Auto' should be a data.frame or a list")
			}
		}
	}

	#========================================
	### Set parameters for Tr and meansParamX
	#========================================
	if(!is.null(X)){
		if(!is.null(Tr)){
			if(is.null(paramTr)){
				if(is.null(meansParamX)){
					paramTr<-matrix(rnorm(nTr*nparamX),nrow=nparamX,ncol=nTr)
					rownames(paramTr)<-paste("p",1:nparamX,sep="")
					colnames(paramTr)<-paste("t",1:nTr,sep="")

					### Community Mean (projector of traits and parameters)
					meansParamX<-matrix(NA,ncol=nparamX,nrow=nsp)
					for(i in 1:nsp){
						meansParamX[i,]<-tcrossprod(Tr[,i],paramTr)
					}
					rownames(meansParamX)<-paste("sp",1:nsp,sep="")
					colnames(meansParamX)<-paste("p",1:nparamX,sep="")
				}else{
					paramTrAll<-array(dim=c(nparamX,nTr,nsp))
					for(i in 1:nsp){
						paramTrAll[,,i]<-meansParamX[i,]%*%matrix(Tr[,i],nrow=1) # matrix() is to account for the case where there is only 1 trait
					}
					### Approximation of paramTr
					paramTr<-apply(paramTrAll,1:2,mean)
				}
			}else{
				if(is.null(meansParamX)){
					### Community Mean (projector of traits and parameters)
					meansParamX<-matrix(NA,ncol=nparamX,nrow=nsp)
					for(i in 1:nsp){
						meansParamX[i,]<-tcrossprod(Tr[,i],paramTr)
					}
					rownames(meansParamX)<-paste("sp",1:nsp,sep="")
					colnames(meansParamX)<-paste("p",1:nparamX,sep="")
				}else{
					### Check if meansParamX and meansParamTr (means calculated using Tr and paramTr) match
					meansParamTr<-matrix(NA,ncol=nparamX,nrow=nsp)
					for(i in 1:nsp){
						meansParamTr[i,]<-tcrossprod(Tr[,i],paramTr)
					}

					if(abs(sum(meansParamX-meansParamTr)) > 10^(-4)){
						warnings("There might be a mismatch between 'meansParamX' and the 'means' calculated using 'Tr' and 'paramTr'")
						warning("The values in 'meansParamX' will be used for the remaining data simulation")
					}

					rownames(meansParamX)<-paste("sp",1:nsp,sep="")
					colnames(meansParamX)<-paste("p",1:nparamX,sep="")
				}
			}
		}else{
			meansParamX<-matrix(rnorm(nparamX),ncol=1)
		}
	}

	#=========================
	### Define outlier species
	#=========================
	outlierWeightVec<-rep(1,nsp)
	if(!is.null(outlierSp)){
		noutlierSp<-round(nsp*outlierSp)
		outlierWeightVec[1:noutlierSp]<-outlierWeight
	}

	#==========================
	### Set parameter for Phylo
	#==========================
	if(!is.null(Phylo)){
		if(is.null(paramPhylo)){
			paramPhylo<-runif(1,0,1)
		}
	}

	#=======================
	### Set parameters for X
	#=======================
	if(!is.null(X)){
		if(is.null(Phylo)){
			### Mean of the community paramX
			if(is.null(paramX)){
				### Sigma matrix of the community
			if(is.null(varX)){
				varX<-chol2inv(chol(rWishart(1,nparamX+2,diag(nparamX))[,,1]))
			}else{
				if(!isSymmetric(varX)){
					stop("'varX' is not a symmetric matrix")
				}
			}
			if(is.null(colnames(varX))){
				colnames(varX)<-paste("p",1:nparamX,sep="")
			}
			if(is.null(rownames(varX))){
				rownames(varX)<-paste("p",1:nparamX,sep="")
			}

			### paramX for each species
			paramX<-matrix(NA,nrow=nsp,ncol=nparamX)
			if(!is.null(Tr)){
				for(i in 1:nsp){
					paramX[i,]<-mvrnorm(1,mu=meansParamX[i,],Sigma=(1/outlierWeightVec[i])*varX)
				}
			}else{
				for(i in 1:nsp){
					paramX[i,]<-mvrnorm(1,mu=meansParamX,Sigma=(1/outlierWeightVec[i])*varX)
				}
			}
				colnames(paramX)<-paste("p",1:nparamX,sep="")
				rownames(paramX)<-paste("sp",1:nsp,sep="")
			}else{
				if(is.null(Tr)){
					if(is.null(meansParamX)){
						meansParamX<-matrix(colMeans(paramX),ncol=1)
					}
				}
				if(is.null(varX)){
					varX<-crossprod(sweep(paramX,2,colMeans(paramX),FUN="-"))
				}
			}
		}else{
			### Sigma matrix of the community
			if(is.null(varX)){
				varX<-chol2inv(chol(rWishart(1,nparamX+2,diag(nparamX))[,,1]))
			}else{
				if(!isSymmetric(varX)){
					stop("'varX' is not a symmetric matrix")
				}
			}
			if(is.null(colnames(varX))){
				colnames(varX)<-paste("p",1:nparamX,sep="")
			}
			if(is.null(rownames(varX))){
				rownames(varX)<-paste("p",1:nparamX,sep="")
			}

			### Weighted phylogeny
			if(paramPhylo>=0){
				wPhylo<-paramPhylo * Phylo + (1-paramPhylo)*diag(nsp)
			}else{
				iPhylo<-cov2cor(chol2inv(chol(Phylo)))
				wPhylo<-(-paramPhylo) * iPhylo + (1+paramPhylo)*diag(nsp)
			}

			varXPhylo <- kronecker(varX,wPhylo);

			if(!is.null(Tr)){
				paramX<-mvrnorm(1,as.vector(meansParamX),varXPhylo)
				paramX<-matrix(paramX,nrow=nsp)
			}else{
				paramX<-mvrnorm(1,rep(meansParamX,nsp),varXPhylo)
				paramX<-matrix(paramX,nrow=nsp,byrow=TRUE)
			}
		}
		#===================
		### Model estimation
		#===================
		EstModel<-tcrossprod(X,paramX)
	}

	#==================================================
	### Set parameters for the latent part of the model
	#==================================================
	if(!is.null(Random)){
		if(is.null(latent)){
			### Object storing the latent variables and the parameters
			latent<-vector("list",length=nRandom)
			dim(latent)<-c(1,ncol(Random))
			colnames(latent)<-paste("random",1:nRandom,sep="")
		}

		if(is.null(paramLatent)){
			paramLatent<-vector("list",length=nRandom)
			dim(paramLatent)<-c(1,ncol(Random))
			names(paramLatent)<-paste("random",1:nRandom,sep="")
		}

		if(is.null(shrinkLocal)){
			### Object storing shrinkLocal
			shrinkLocal<-vector("list",length=nRandom)
			dim(shrinkLocal)<-c(1,ncol(Random))
			colnames(shrinkLocal)<-paste("random",1:nRandom,sep="")
		}else{
			if(length(shrinkLocal)!=nRandom){
				stop("'shrinkLocal' needs to have the same length Random")
			}
		}

		if(is.null(paramShrinkGlobal)){
			### Object storing shrinkLocal
			paramShrinkGlobal<-vector("list",length=nRandom)
			dim(paramShrinkGlobal)<-c(1,ncol(Random))
			colnames(paramShrinkGlobal)<-paste("random",1:nRandom,sep="")
		}else{
			if(length(paramShrinkGlobal)!=nRandom){
				stop("'paramShrinkGlobal' needs to have the same length Random")
			}
		}

		latentCov<-array(dim=c(nsp,nsp,nRandom))
		dimnames(latentCov)[[1]]<-paste("sp",1:nsp,sep="")
		dimnames(latentCov)[[2]]<-paste("sp",1:nsp,sep="")
		dimnames(latentCov)[[3]]<-paste("random",1:nRandom,sep="")

		### Number of levels in each random effect considered
		nLevelRandom<-mapply(nlevels,Random)
		### Construct the model
		for(i in 1:nRandom){
			if(!is.null(latent[[1,i]]) & !is.null(paramLatent[[1,i]])){
				nLatent <- ncol(latent[[1,i]])

				### Check if the levels match between Random and latent
				matchLev <- which(levels(Random[,i]) %in% rownames(latent[[1,i]]))
				locLev <- which(rownames(latent[[1,i]]) %in% Random[[matchLev,i])

				### Build new latent variables
				latentNew<-matrix(NA,nrow=nLevelRandom[i],ncol=nLatent)
				### If none of the levels match
				if(length(matchLev) == 0){
					latentNew <- matrix(rnorm(nLevelRandom[i]*nLatent),nrow=nLevelRandom[i],ncol=nLatent)
				}

				### If the levels match
				if(length(matchLev) == nLevelRandom){
					latentNew <- latent[[1,i]][locLev,]
				}

				### If only some of the levels match
				if(length(matchLev) > 0 && length(matchLev) < nLevelRandom[i]){
					latentNew[matchLev,] <- latent[[1,i]][locLev,]
					latentNew[-matchLev,] <- matrix(norm((nLevelRandom[i]-length(matchLev))*nLatent),nrow=nLevelRandom[i],ncol=nLatent)
				}
				### Assign latentNew to latent
				latent[[1, i]] <- latentNew[[1, i]]
			}

			### Define the latent variables and their associated parameters
			if(is.null(latent[[1,i]])){
				if(is.null(paramLatent[[1,i]])){
					if(is.null(paramShrinkGlobal[[1,i]])){
						if(is.null(shrinkLocal[[1,i]])){
							### Define the number of latent variables in the model (it will vary with the number of random effects)
							nLatent<-sample(2:5,nRandom)
						}else{
							nLatent<-sapply(shrinkLocal,ncol)
						}
					}else{
						nLatent<-sapply(paramShrinkGlobal,length)
					}
				}else{
					if(is.null(paramShrinkGlobal[[1,i]])){
						if(is.null(shrinkLocal[[1,i]])){
							### Define the number of latent variables in the model (it will vary with the number of random effects)
							nLatent<-sapply(paramLatent[[1,i]],ncol)
						}else{
							nLatent<-sapply(shrinkLocal,ncol)
						}
					}else{
						nLatent<-sapply(paramShrinkGlobal,length)
					}
				}
				### Construct latent variables
				latent[[1,i]]<-matrix(rnorm(nLevelRandom[i]*nLatent[i]),nrow=nLevelRandom[i],ncol=nLatent[i])
				rownames(latent[[1,i]])<-paste("lev",1:nLevelRandom[i],sep="")
				colnames(latent[[1,i]])<-paste("pl",1:nLatent[i],sep="")
			}

			### Define paramLatent and their associated parameters
			if(is.null(paramLatent[[1,i]])){
				if(is.null(shrinkLocal[[1,i]])){
					shrinkLocal[[1,i]]<-matrix(rgamma(nsp*nLatent[i],1.5,1.5),nsp,nLatent[i])
				}else{
					nLatentTest<-sapply(shrinkLocal,ncol)
					if(!all(nLatentTest%in%nLatent)){
						stop("'shrinkLocal' needs to have the same number levels as Random")
					}
				}

				if(is.null(paramShrinkGlobal[[1,i]])){
					paramShrinkGlobal[[1,i]]<-c(rgamma(1,50,1),rgamma(nLatent[i]-1,50,1))
				}else{
					nLatentTest<-sapply(paramShrinkGlobal,length)
					if(!all(nLatentTest%in%nLatent)){
						stop("'paramShrinkGlobal' needs to have the same number levels as Random")
					}
				}

				### Construct latent variables
				shrinkGlobal<-cumprod(paramShrinkGlobal[[1,i]])
				normSD<-as.vector(1/matrix(rep(shrinkGlobal,nsp),nrow=nsp,byrow=TRUE)*shrinkLocal[[1,i]])

				paramLatent[[1,i]]<-matrix(rnorm(length(normSD),0,normSD),nsp,nLatent[i])
			}
			if(!is.null(X)){
				### Add the random effect to the estimated model
				EstModel<-EstModel+tcrossprod(latent[[1,i]][Random[,i],],paramLatent[[1,i]])
			}else{
				### Estimated model
				EstModel<-tcrossprod(latent[[1,i]][Random[,i],],paramLatent[[1,i]])
			}

			### Calculate the variance-covariance matrix from the latent variables
			latentCov[,,i]<-tcrossprod(paramLatent[[1,i]])
		}
	}

	if(!is.null(X)){
		if(!is.null(spCor)){
			EstModel<-EstModel+mvrnorm(n=nsite,mu=rep(0,nsp),Sigma=spCor)
		}
	}else{
		if(!is.null(spCor)){
			EstModel<-mvrnorm(n=nsite,mu=rep(0,nsp),Sigma=spCor)
		}
	}

	#===========================================================
	### Set parameters for the autocorrelation part of the model
	#===========================================================
	if(!is.null(Auto)){
		#__________________________________________________________
		### Calculate the distance between sampled within one level
		#__________________________________________________________
		### Number of levels in each auto effect considered
		nLevelAuto<-sapply(Auto,function(x) nlevels(x[,1]))
		AutoDist<-vector("list",length=nAuto)

		for(i in 1:nAuto){
			nAutoCoord<-ncol(Auto[[i]])-1
			AutoCoordMean<-matrix(NA,nrow=nLevelAuto[i],ncol=nAutoCoord)

			for(j in 1:nAutoCoord){
				AutoCoordMean[,j]<-tapply(Auto[[i]][,j+1],Auto[[i]][,1],mean)
			}
			AutoDist[[i]]<-dist(AutoCoordMean)
		}

		if(is.null(latentAuto)){
			### Construct latentAuto (no values yet)
			latentAuto<-vector("list",length=nAuto)
			dim(latentAuto)<-c(1,nAuto)
			colnames(latentAuto)<-paste("auto",1:nAuto,sep="")
		}

		if(is.null(paramLatentAuto)){
			paramLatentAuto<-vector("list",length=nAuto)
			dim(paramLatentAuto)<-c(1,nAuto)
			names(paramLatentAuto)<-paste("auto",1:nAuto,sep="")
		}

		if(is.null(shrinkLocalAuto)){
			### Object storing shrinkLocalAuto
			shrinkLocalAuto<-vector("list",length=nAuto)
			dim(shrinkLocalAuto)<-c(1,nAuto)
			colnames(shrinkLocalAuto)<-paste("auto",1:nAuto,sep="")
		}else{
			if(length(shrinkLocalAuto)!=nAuto){
				stop("'shrinkLocalAuto' needs to have the same length Auto")
			}
		}

		if(is.null(paramShrinkGlobalAuto)){
			### Object storing paramShrinkGlobalAuto
			paramShrinkGlobalAuto<-vector("list",length=nAuto)
			dim(paramShrinkGlobalAuto)<-c(1,nAuto)
			colnames(paramShrinkGlobalAuto)<-paste("auto",1:nAuto,sep="")
		}else{
			if(length(paramShrinkGlobalAuto)!=nAuto){
				stop("'paramShrinkGlobalAuto' needs to have the same length Auto")
			}
		}

		### Number of levels in each autocorrelated random effect considered
		nLevelAuto<-mapply(nlevels,Auto)

		### Construct the model
		for(i in 1:nAuto){
			if(!is.null(latentAuto[[1,i]]) & !is.null(paramLatentAuto[[1,i]])){
				nLatent <- ncol(latentAuto[[1,i]])

				### Check if the levels match between Auto and latentAuto
				matchLev <- which(levels(Auto[,i]) %in% rownames(latentAuto[[1,i]]))
				locLev <- which(rownames(latentAuto[[1,i]]) %in% Auto[[matchLev,i])

				### Build new latentAuto variables
				latentNew<-matrix(NA,nrow=nLevelAuto[i],ncol=nLatent)
				### If none of the levels match
				if(length(matchLev) == 0){
					latentNew <- matrix(rnorm(nLevelAuto[i]*nLatent),nrow=nLevelAuto[i],ncol=nLatent)
				}

				### If the levels match
				if(length(matchLev) == nLevelAuto){
					latentNew <- latentAuto[[1,i]][locLev,]
				}

				### If only some of the levels match
				if(length(matchLev) > 0 && length(matchLev) < nLevelAuto[i]){
					latentNew[matchLev,] <- latentAuto[[1,i]][locLev,]
					latentNew[-matchLev,] <- matrix(norm((nLevelAuto[i]-length(matchLev))*nLatent),nrow=nLevelAuto[i],ncol=nLatent)
				}
				### Assign latentNew to latentAuto
				latentAuto[[1, i]] <- latentNew[[1, i]]
			}

			### Define the latent variables and their associated parameters
			if(is.null(latentAuto[[1,i]])){
				if(is.null(paramLatentAuto[[1,i]])){
					if(is.null(paramShrinkGlobalAuto[[1,i]])){
						if(is.null(shrinkLocalAuto[[1,i]])){
							### Define the number of autocorrelated latent variables in the model (it will vary with the number of random effects)
							nLatentAuto<-sample(2:5,nAuto)
						}else{
							nLatentAuto<-sapply(shrinkLocalAuto,ncol)
						}
					}else{
						nLatentAuto<-sapply(paramShrinkGlobalAuto,length)
					}
				}else{
					if(is.null(paramShrinkGlobalAuto[[1,i]])){
						if(is.null(shrinkLocalAuto[[1,i]])){
							### Define the number of autocorrelated latent variables in the model (it will vary with the number of random effects)
							nLatentAuto<-sapply(paramLatentAuto[[1,i]],ncol)
						}else{
							nLatentAuto<-sapply(shrinkLocalAuto,ncol)
						}
					}else{
						nLatentAuto<-sapply(paramShrinkGlobalAuto,length)
					}
				}

				#--------------------------------------
				### Define paramAuto if it is not given
				#--------------------------------------
				if(is.null(paramAuto)){
					paramAuto<-vector("list",length=nAuto)
					for(i in 1:nAuto){
						paramAuto[[i]]<-runif(nLatentAuto[i])*max(AutoDist[[1]])
					}
				}else{
					if(length(paramAuto)!=nAuto){
						stop("The length of 'paramAuto' should be the same as the length of Auto")
					}
					if(all(sapply(paramAuto,length)!=nLatentAuto)){
						stop("Each part of 'paramAuto' should have the same length as the number of latent variables in each part of Auto")
					}
				}
				#-----------------------------------------------------------
				### Construct object to weighted the distance with paramAuto
				#-----------------------------------------------------------
				### Use Exponential function
				wAutoDist<-vector("list",length=nAuto)
				for(i in 1:nAuto){
					wAutoDist[[i]]<-vector("list",length=nLatentAuto[i])
				}

				for(i in 1:nAuto){
					AutoDistMat<-as.matrix(AutoDist[[i]])
					for(j in 1:nLatentAuto[i]){
						wAutoDist[[i]][[j]]<-exp(-AutoDistMat/paramAuto[[i]][j])
					}
				}

				#__________________________________
				### Autocorrelated latent variables
				#__________________________________
				for(i in 1:nAuto){
					latentAuto[[1,i]]<-matrix(NA,nrow=nLevelAuto[i],ncol=nLatentAuto[i])
					for(j in 1:nLatentAuto[i]){
						latentAuto[[1,i]][,j]<-rmvnorm(1,rep(0,nrow(wAutoDist[[i]][[j]])),wAutoDist[[i]][[j]])
					}
					rownames(latentAuto[[1,i]])<-levels(Auto[[i]][,1])
					colnames(latentAuto[[1,i]])<-paste("latentAuto",1:ncol(latentAuto[[1,i]]),sep="")
				}
			}

			### Define paramLatentAuto and their associated parameters
			if(is.null(paramLatentAuto[[1,i]])){
				if(is.null(shrinkLocalAuto[[1,i]])){
					shrinkLocalAuto[[1,i]]<-matrix(rgamma(nsp*nLatentAuto[i],1.5,1.5),nsp,nLatentAuto[i])
				}else{
					nLatentAutoTest<-sapply(shrinkLocalAuto,ncol)
					if(!all(nLatentAutoTest%in%nLatentAuto)){
						stop("'shrinkLocalAuto' needs to have the same number levels as Auto")
					}
				}

				if(is.null(paramShrinkGlobalAuto[[1,i]])){
					paramShrinkGlobalAuto[[1,i]]<-c(rgamma(1,50,1),rgamma(nLatentAuto[i]-1,50,1))
				}else{
					nLatentAutoTest<-sapply(paramShrinkGlobalAuto,length)
					if(!all(nLatentAutoTest%in%nLatentAuto)){
						stop("'paramShrinkGlobalAuto' needs to have the same number levels as Random")
					}
				}

				### Construct latent variables
				shrinkGlobalAuto<-cumprod(paramShrinkGlobalAuto[[1,i]])
				normSD<-as.vector(1/matrix(rep(shrinkGlobalAuto,nsp),nrow=nsp,byrow=TRUE)*shrinkLocalAuto[[1,i]])

				paramLatentAuto[[1,i]]<-matrix(rnorm(length(normSD),0,normSD),nsp,nLatentAuto[i])
			}

			if(!is.null(Random)){
				if(!is.null(X)){
					### Add the auto effect to the estimated model
					EstModel<-EstModel+tcrossprod(latentAuto[[1,i]][Auto[[i]][,1],],paramLatentAuto[[1,i]])
				}
			}else{
				if(!is.null(X)){
					### Add the auto effect to the estimated model
					EstModel<-EstModel+tcrossprod(latentAuto[[1,i]][Auto[[i]][,1],],paramLatentAuto[[1,i]])
				}else{
					### Estimated model
					EstModel<-tcrossprod(latentAuto[[1,i]][Auto[[i]][,1],],paramLatentAuto[[1,i]])
				}
			}
		}
	}

	#======================================
	### Construct species occurrence matrix
	#======================================

	if(family=="gaussian"){
		Y<-matrix(rnorm(nsite*nsp,mean=as.vector(EstModel),sd=paramDist),nrow=nsite,ncol=nsp)
	}
	if(family=="nbinomial"){
		Y<-matrix(rnbinom(nsite*nsp,mu=as.vector(EstModel),size=paramDist),nrow=nsite,ncol=nsp)
	}
	if(family=="probit"){
		Ylatent<-matrix(rnorm(nsite*nsp,mean=as.vector(EstModel),sd=1),nrow=nsite,ncol=nsp)
		Y<-Ylatent
		Y[Y>0]<-1 # * much faster than ifelse()
		Y[Y<0]<-0
	}
	if(family=="poisson"){
		Y<-matrix(rpois(nsite*nsp,lambda=as.vector(exp(EstModel))),nrow=nsite,ncol=nsp)
	}

	rownames(Y)<-paste("site",1:nsite,sep="")
	colnames(Y)<-paste("sp",1:nsp,sep="")

	#=============================
	### Return the results objects
	#=============================
	### Data
	if(is.null(Auto)){
		if(!is.null(X)){
			if(is.null(Phylo)){
				if(!is.null(Random)){
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Random=Random)
						attributes(data)<-list(names=c("Y","X","Tr","Random"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Random=Random)
						attributes(data)<-list(names=c("Y","X","Random"),Ypattern="sp")
					}
				}else{
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr))
						attributes(data)<-list(names=c("Y","X","Tr"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X))
						attributes(data)<-list(names=c("Y","X"),Ypattern="sp")
					}
				}
			}else{
				if(!is.null(Random)){
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo),Random=Random)
						attributes(data)<-list(names=c("Y","X","Tr","Phylo","Random"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo),Random=Random)
						attributes(data)<-list(names=c("Y","X","Phylo","Random"),Ypattern="sp")
					}
				}else{
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo))
						attributes(data)<-list(names=c("Y","X","Tr","Phylo"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo))
						attributes(data)<-list(names=c("Y","X","Phylo"),Ypattern="sp")
					}
				}
			}
		}else{
			if(!is.null(Random)){
				data<-list(Y=as.matrix(Y),Random=Random)
				attributes(data)<-list(names=c("Y","Random"),Ypattern="sp")
			}
		}
	}else{
		if(!is.null(X)){
			if(is.null(Phylo)){
				if(!is.null(Random)){
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Random=Random,Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Tr","Random","Auto"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Random=Random,Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Random","Auto"),Ypattern="sp")
					}
				}else{
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Tr","Auto"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Auto"),Ypattern="sp")
					}
				}
			}else{
				if(!is.null(Random)){
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo),Random=Random,Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Tr","Phylo","Random","Auto"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo),Random=Random,Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Phylo","Random","Auto"),Ypattern="sp")
					}
				}else{
					if(!is.null(Tr)){
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Tr=as.matrix(Tr),Phylo=as.matrix(Phylo),Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Tr","Phylo","Auto"),Ypattern="sp")
					}else{
						data<-list(Y=as.matrix(Y),X=as.matrix(X),Phylo=as.matrix(Phylo),Auto=Auto)
						attributes(data)<-list(names=c("Y","X","Phylo","Auto"),Ypattern="sp")
					}
				}
			}
		}else{
			if(!is.null(Random)){
				data<-list(Y=as.matrix(Y),Random=Random,Auto=Auto)
				attributes(data)<-list(names=c("Y","Random","Auto"),Ypattern="sp")
			}else{
				data<-list(Y=as.matrix(Y),Auto=Auto)
				attributes(data)<-list(names=c("Y","Auto"),Ypattern="sp")
			}
		}

	}
	class(data)<-"HMSCdata"

	### Parameters
	if(is.null(Auto)){
		if(!is.null(X)){
			if(is.null(Phylo)){
				if(!is.null(Random)){
					if(!is.null(Tr)){
						allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
					}else{
						allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
					}
				}else{
					if(!is.null(spCor)){
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor)
						}
					}else{
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec)
						}
					}
				}
			}else{
				if(!is.null(Random)){
					if(!is.null(Tr)){
						allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
					}else{
						allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
					}
				}else{
					if(!is.null(spCor)){
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor)
						}
					}else{
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec)
						}
					}
				}
			}
		}else{
			if(!is.null(Random)){
				allparam<-list(outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal)
			}else{
				if(!is.null(spCor)){
					allparam<-list(outlier=outlierWeightVec,spCor=spCor)
				}
			}
		}
	}else{
		if(!is.null(X)){
			if(is.null(Phylo)){
				if(!is.null(Random)){
					if(!is.null(Tr)){
						allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
					}else{
						allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
					}
				}else{
					if(!is.null(spCor)){
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}
					}else{
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}
					}
				}
			}else{
				if(!is.null(Random)){
					if(!is.null(Tr)){
						allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
					}else{
						allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
					}
				}else{
					if(!is.null(spCor)){
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}
					}else{
						if(!is.null(Tr)){
							allparam<-list(paramX=paramX,paramTr=paramTr, meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}else{
							allparam<-list(paramX=paramX,meansParamX=meansParamX,varX=varX,precX=solve(varX),paramPhylo=paramPhylo,outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
						}
					}
				}
			}
		}else{
			if(!is.null(Random)){
				allparam<-list(outlier=outlierWeightVec,latent=latent,paramLatent=paramLatent,latentCov=latentCov,shrinkLocal=shrinkLocal,paramShrinkGlobal=paramShrinkGlobal,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
			}else{
				if(!is.null(spCor)){
					allparam<-list(outlier=outlierWeightVec,spCor=spCor,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
				}else{
					allparam<-list(outlier=outlierWeightVec,latentAuto=latentAuto,paramLatentAuto=paramLatentAuto,paramAuto=paramAuto,shrinkLocalAuto=shrinkLocalAuto,paramShrinkGlobalAuto=paramShrinkGlobalAuto)
				}
			}
		}
	}

	class(allparam)<-"HMSCparam"

	### Make sure allparam are all matrices
	allparam<-lapply(allparam,as.matrix)

	### Final result object
	if(family=="gaussian"){
		res<-list(data=data,param=allparam,sd=paramDist,probMat=EstModel)
	}
	if(family=="nbinomial"){
		res<-list(data=data,param=allparam,size=paramDist,probMat=EstModel)
	}
	if(family=="probit"){
		res<-list(data=data,param=allparam,probMat=pnorm(EstModel))
	}
	if(family=="poisson"){
		res<-list(data=data,param=allparam,probMat=exp(EstModel))
	}
	class(res)<-"communitySimul"

	return(res)
}
