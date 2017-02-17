#' @title Name result
#'
#' @description Name all parts of the result object. This is meant to be an internal function.
#'
#' @param data An object of the class \code{HMSCdata}.
#' @param priors An object of the class \code{HMSCprior}.
#' @param result An hmsc object without names on any part of the object.
#' @param listLev A numeric value defining on which level of the list names should be added.
#' @param niter An integer defining the number of iterations used to perform the analysis.
#' @param nburn An integer defining the number of burning iterations used to perform the analysis.
#' @param thin An integer defining thinning in the analysis.
#' 
#' @return
#'
#' The result object with name for all dimensions and parts of lists
#'
#' @author F. Guillaume Blanchet
#'
#' @keywords utilities
#' @export
nameResult <-
function(data,priors,result,niter,nburn,thin){
	
	### Number of parameters estimated
	nParam<-length(names(result[[1]]))
	
	if(attr(priors,"distr")=="probit"){
		if(nParam==2){
			#==============
			### only Random
			#==============
			if(all(names(result[[1]])==c("paramLatent","latent"))){
				result<-nameRandom(data,result,niter,nburn,thin,listLev=1)
			}
		}
		if(nParam==3){
			#=========
			### only X
			#=========
			if(all(names(result[[1]])==c("paramX","meansParamX","varX"))){
				result<-nameX(data,result,niter,nburn,thin)
			}
			
			#===========
			### X and Tr
			#===========
			if(all(names(result[[1]])==c("paramX","paramTr","varX"))){
				result<-nameXTr(data,result,niter,nburn,thin)
			}
			
			#============
			### only Auto
			#============
			if(all(names(result[[1]])==c("paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameAuto(data,result,niter,nburn,thin,listLev=1)
			}
		}
		
		if(nParam==4){
			#==============
			### X and Phylo
			#==============
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","paramPhylo"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=4)
			}
		
			#==================
			### X, Tr and Phylo
			#==================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","paramPhylo"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=4)
			}
		}

		if(nParam==5){
			#===============
			### X and Random
			#===============
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","paramLatent","latent"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=4)
			}
		
			#===================
			### X, Tr and Random
			#===================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","paramLatent","latent"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=4)
			}
		
			#==================
			### Random and Auto
			#==================
			if(all(names(result[[1]])==c("paramLatent","latent","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameRandom(data,result,niter,nburn,thin,listLev=1)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=3)
			}
		}
		
		if(nParam==6){
			#=============
			### X and Auto
			#=============
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=4)
			}
			
			#=================
			### X, Tr and Auto
			#=================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=4)
			}
			
			#======================
			### X, Phylo and Random
			#======================
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","paramPhylo","paramLatent","latent"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=4)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=5)
			}
		
			#==========================
			### X, Tr, Phylo and Random
			#==========================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","paramPhylo","paramLatent","latent"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=4)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=5)
			}
		}
		
		if(nParam==7){
			#====================
			### X, Phylo and Auto
			#====================
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","paramPhylo","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=4)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=5)
			}
			
			#========================
			### X, Tr, Phylo and Auto
			#========================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","paramPhylo","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=4)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=5)
			}
		}
		
		if(nParam==8){
			#=====================
			### X, Random and Auto
			#=====================
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","paramLatent","latent","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=4)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=6)
			}
		}
	}
	
	### For Gaussian data
	if(attr(priors,"distr")=="gaussian"){
		if(nParam==3){
			#==============
			### only Random
			#==============
			if(all(names(result[[1]])==c("varNormal","paramLatent","latent"))){
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=1)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=2)
			}
		}
		if(nParam==4){
			#=========
			### only X
			#=========
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
			}
			
			#===========
			### X and Tr
			#===========
			if(all(names(result[[1]])==c("paramX","paramTr","varX","varNormal"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
			}
			
			#============
			### only Auto
			#============
			if(all(names(result[[1]])==c("varNormal","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=1)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=2)
			}
		}
		
		if(nParam==5){
			#==============
			### X and Phylo
			#==============
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal","paramPhylo"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=5)
			}
		
			#==================
			### X, Tr and Phylo
			#==================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","varNormal","paramPhylo"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=5)
			}
		}

		if(nParam==6){
			#===============
			### X and Random
			#===============
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal","paramLatent","latent"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=5)
			}
		
			#===================
			### X, Tr and Random
			#===================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","varNormal","paramLatent","latent"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=5)
			}
			
			#==================
			### Random and Auto
			#==================
			if(all(names(result[[1]])==c("varNormal","paramLatent","latent","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=1)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=2)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=4)
			}
		}
		
		if(nParam==7){
			#=============
			### X and Auto
			#=============
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=5)
			}
			
			#=================
			### X, Tr and Auto
			#=================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","varNormal","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=5)
			}
			
			#======================
			### X, Phylo and Random
			#======================
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal","paramPhylo","paramLatent","latent"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=5)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=6)
			}
		
			#==========================
			### X, Tr, Phylo and Random
			#==========================
			if(all(names(result[[1]])==c("paramX","paramTr","varX","varNormal","paramPhylo","paramLatent","latent"))){
				result<-nameXTr(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=5)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=6)
			}
		}
		
		if(nParam==8){
			#====================
			### X, Phylo and Auto
			#====================
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal","paramPhylo","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-namePhylo(data,result,niter,nburn,thin,listLev=5)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=6)
			}
		}
		
		if(nParam==9){
			#=====================
			### X, Random and Auto
			#=====================
			if(all(names(result[[1]])==c("paramX","meansParamX","varX","varNormal","paramLatent","latent","paramLatentAuto","latentAuto","paramAuto"))){
				result<-nameX(data,result,niter,nburn,thin)
				result<-nameResidVar(data,result,niter,nburn,thin,listLev=4)
				result<-nameRandom(data,result,niter,nburn,thin,listLev=5)
				result<-nameAuto(data,result,niter,nburn,thin,listLev=7)
			}
		}
	}
	
	### Name the result
	if(length(result)==2){
		names(result)[[1]]<-"burning"
		names(result)[[2]]<-"estimation"
	}
	if(length(result)==1){
		names(result)[[1]]<-"estimation"
	}
	
	return(result)
}
