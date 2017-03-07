#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcPoissonX.h"

using namespace arma ;
using namespace Rcpp ;

//' @title Markov chain Monte Carlo algorithms.
//'
//' @description These functions are meant to be used internally or by advanced users as they generate raw (somewhat unnamed output)
//'
//' @param Y Raw site by species community matrix.
//' @param Ylatent Model site by species community matrix after the link function is applied.
//' @param X Site by covariates matrix.
//' @param Tr Trait by species matrix.
//' @param Phylo Square correlation (with values ranging from -1 to 1) matrix with as many rows (and columns) as there are species in \code{Y}.
//' @param iPhylo Square correlation matrix with as many rows (and columns) as there are species in \code{Y}. Inverse of \code{Phylo}.
//' @param Auto List with as many levels as there are autocorrelated factors in the analysis. Each level of the list contains a set of coordinates associated to each level of each factor.
//' @param RandomAuto A data.frame with all the autocorrelated factors (organized as columns) in the analysis. The number of columns in \code{RandomAuto} equals the number of levels in the \code{Auto}.
//' @param Random A data.frame with all the non-autocorrelated factors (organized as columns) in the analysis.
//' @param paramX Species by covariates matrix of regression parameters.
//' @param meansParamX Matrix with a single column and as many rows as covariates. This matrix defines how an average species react to the covariates.
//' @param paramTr Covariate by traits matrix of regression parameters.
//' @param precX Square (precision) matrix with as many rows (and columns) as there are covariates. The inverse of this matrix is a variance matrix which explains how a typical species vary in its response to one (diagonal) or a pair (off-diagonal) of covariates.
//' @param paramPhylo Numeric. Single parameter describing the importance of phylogeny in the analysis.
//' @param residVar Vector of parameters with as many entries as there are species. Each values in this vector describes how much variation is not explained by the model.
//' @param latent List. Each level of the list includes a matrix of latent (unsampled) variables. There are as many level in the list as there are columns in \code{Random}. The number of rows in each matrix corresponds to the number of levels in each factor in \code{Random}.
//' @param paramLatent List. Each level of the list includes a matrix of parameters for the latent (unsampled) variables. There are as many level in the list as there are columns in \code{Random}. Each matrix has as many rows as there are species in \code{Y} and as many columns as there are latent variables in \code{latent}.
//' @param shrinkLocal List. Each level of the list includes a matrix of values describing how much shrinkage should be imposed on each parameter of the latent (unsampled) variables (\code{paramLatent}). There are as many level in the list as there are columns in \code{Random}. The size of each matrix is the same as the size of \code{paramLatent}.
//' @param paramShrinkGlobal List. Each level of the list includes a vector of values describing how much shrinkage should be imposed on each latent variables (\code{latent}) as a whole. There are as many level in the list as there are columns in \code{Random}. The number of elements in each vector equals the number of variables in \code{latent}.
//' @param paramAuto List.  Each level of the list includes a vector of parameters, one for each autocorrelated latent variables in \code{latentAuto}.  There are as many level in the list as there are columns in \code{Auto}.
//' @param latentAuto List. Each level of the list includes a matrix of autocorrelated latent (unsampled) variables. There are as many level in the list as there are columns in \code{Auto}. The number of rows in each matrix corresponds to the number of levels in each autocorrelated factor in \code{Auto}.
//' @param paramLatentAuto List. Each level of the list includes a matrix of parameters for the autocorrelated latent (unsampled) variables. There are as many level in the list as there are columns in \code{Auto}. Each matrix has as many rows as there are species in \code{Y} and as many columns as there are autocorrelated latent variables in \code{latentAuto}.
//' @param shrinkLocalAuto List. Each level of the list includes a matrix of values describing how much shrinkage should be imposed on each parameter of the autocorrelated latent (unsampled) variables (\code{paramLatentAuto}). There are as many level in the list as there are columns in \code{Auto}. The size of each matrix is the same as the size of \code{paramLatentAuto}.
//' @param paramShrinkGlobalAuto List. Each level of the list includes a vector of values describing how much shrinkage should be imposed on each autocorrelated latent variables (\code{latent}) as a whole. There are as many level in the list as there are columns in \code{Auto}. The number of elements in each vector equals the number of variables in \code{latentAuto}.
//' @param priorMeansParamX Prior. Matrix with a single column and as many rows as covariates. It describes the prior values associated to meansParamX.
//' @param priorVarMeansParamX square matrix of parameters defining how \code{meansParamX} varies.
//' @param priorParamTr Prior. Matrix of priors defining how descriptors (rows) characterizes traits (columns).
//' @param priorVarTr Prior. Symmetric covariance matrix of prior. Each dimension of this matrix is be equal to the number of traits.
//' @param priorVarXScaleMat Square matrix with as many rows (and columns) as there are covariates. Prior. Scale matrix used to sample \code{precX} from an inverse Wishart distribution.
//' @param priorVarXDf Numeric. Prior. Number of degrees of freedoms used to sample \code{precX} from an inverse Wishart distribution.
//' @param priorResidVarScale Numeric. Prior. Scale parameter of a gamma distribution from which \code{residVar} is sampled.
//' @param priorResidVarShape Numeric. Prior. Shape parameter of a gamma distribution from which \code{residVar} is sampled.
//' @param priorParamPhyloWeight Prior. Matrix with a single column defining the importance each potential phylogenetic weight (given in \code{priorParamPhyloGrid}) may have. The size of \code{priorParamPhyloWeight} needs to be equal to the size of \code{priorParamPhyloGrid}.
//' @param priorParamPhyloGrid Prior. Matrix with a single column defining the potential values \code{paramPhylo} can take. The size of \code{priorParamPhyloGrid} needs to be equal to the size of \code{priorParamPhyloWeight}.
//' @param priorParamAutoWeight Prior. Matrix with a single column defining the importance each value in \code{priorParamAutoDist} may have. The size of \code{priorParamAutoWeight} needs to be equal to the size of \code{priorParamAutoGrid}.
//' @param priorParamAutoDist Prior. Matrix with a single column defining the potential values given in \code{paramAuto} can take. The size of \code{priorParamAutoGrid} needs to be equal to the size of \code{priorParamAutoWeight}.
//' @param priorShrinkLocal Numeric. Hyperparameter of a gamma distribution defining the local shrinkage for latent variables.
//' @param priorShrinkOverallShape Numeric. Shape of a gamma distribution.
//' @param priorShrinkOverallScale Numeric. Scale of a gamma distribution.
//' @param priorShrinkSpeedShape Numeric. Shape of a gamma distribution.
//' @param priorShrinkSpeedScale Numeric. Scale of a gamma distribution.
//' @param adapt Vector of two parameters defining whether the number of latent and autocorrelated latent variables should be modified. The function used is : 1/exp(adapt[1]+(adapt[2]*niter)).
//' @param redund Vector of two parameters defining (1) the proportion of redundant elements within latent and autocorrelated latent variables and (2) the amount of error that can be accounted for given the proportion of redundancy.
//' @param nAuto Numeric. Number of autocorrelated random effect.
//' @param nAutoLev Vector defining the number of levels for each autocorrelated random effect. The length of this vector equals \code{nAuto}.
//' @param nLatentAuto Vector defining the number of autocorrelated latent variables for each autocorrelated random effect. The length of this vector equals \code{nAuto}.
//' @param nRandom Numeric. Number of non-autocorrelated random effect.
//' @param nRandomLev Vector defining the number of levels for each non-autocorrelated random effect. The length of this vector equals \code{nRandom}.
//' @param nLatent Vector defining the number of non-autocorrelated latent variables for each non-autocorrelated random effect. The length of this vector equals \code{nRandom}.
//' @param nsp Numeric. Number of species in \code{Y}.
//' @param nsite Numeric. Number of site sampled.
//' @param nparamX Numeric. Number of covariates.
//' @param nTr Numeric. Number of traits.
//' @param nparamPhylo Numeric. Number of values in \code{priorParamPhyloGrid}.
//' @param npriorParamAuto Numeric. Number of values in \code{priorParamAutoGrid}.
//' @param niter Numeric. Number of iterations performed.
//' @param nburn Numeric. Of the number of iterations performed, number of discarded iterations. Burning phase of the model estimations.
//' @param thin Numeric. Thinning. Iterations saved every time a multiple of the given value is an integer.
//' @param verbose Logical or numeric. If \code{TRUE}, the number of iterations are printed on the screen four times. If \code{FALSE}, nothing is printed on the screen. If a positive integer is given, the number of iterations is printed on the screen every time a multiple of the given value is an integer.
//'
//' @export
//[[Rcpp::export]]
RcppExport SEXP mcmcPoissonX(arma::mat& Y,
							arma::mat& Ylatent,
							arma::mat& X,
							arma::mat& paramX,
							arma::mat& meansParamX,
							arma::mat& precX,
							arma::vec& residVar,
							arma::mat& priorMeansParamX,
							arma::mat& priorVarXScaleMat,
							double priorVarXDf,
							double priorResidVarScale,
							double priorResidVarShape,
							double nsp,
							int nsite,
							int nparamX,
							int niter,
							int nburn,
							int thin,
							int verbose){

	// Define various objects
	mat EstModel(nsite,nsp);
	mat Yresid(nsite,nsp);
	mat meansParamXRep(nsp, nparamX);

	// Define the result objects for burning
	mat meansParamXBurn(nparamX,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	mat meansParamXEst(nparamX,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);

	// Define a burning counters
	int countBurn;
	countBurn = 0;

	// Define an estimation counters
	int countEst;
	countEst = 0;

	// Gibbs sampler
	for (int i = 0; i < niter; i++) {
		// Calculate the model estimate

		EstModel = X*trans(paramX);
//Rprintf("1\n");
		// Sample Y latent
		Ylatent = sampleYlatentPoisson(Y, Ylatent, EstModel, residVar, nsp, nsite);
//Rprintf("2\n");

		// Update paramX
		meansParamXRep = trans(repmat(meansParamX,1,nsp));
		paramX = updateParamX(Ylatent,X,meansParamXRep,precX, paramX, nsp, nsite, nparamX);

		// Update precX
		precX = updatePrecX(meansParamX,priorVarXScaleMat, priorVarXDf, paramX, precX, nsp);

		// Update meanparamX
//		Rprintf("4.5\n");
		meansParamX = updateMeansParamX(priorMeansParamX, priorVarXScaleMat, priorVarXDf, paramX, meansParamX, precX, nsp, nparamX);
//		Rprintf("5\n");

		if(i<nburn && i%thin==0){
			// Save burning results
			meansParamXBurn.col(countBurn) = meansParamX;
			varXBurn.slice(countBurn) = precX.i();
			paramXBurn.slice(countBurn) = paramX;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			meansParamXEst.col(countEst) = meansParamX;
			varXEst.slice(countEst) = precX.i();
			paramXEst.slice(countEst) = paramX;

			// Counter
			countEst++;
		}

		//Print status of MCMC run
		if (i>1 && i%verbose==0) {
			Rprintf("iteration %d\n",i);
		}

	}

	// Return a list of results
	return Rcpp::List::create(
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXBurn)),
								   Rcpp::Named("meansParamX", wrap(trans(meansParamXBurn))),
							 	   Rcpp::Named("varX", wrap(varXBurn))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("meansParamX", wrap(trans(meansParamXEst))),
							 	   Rcpp::Named("varX", wrap(varXEst))));

}
