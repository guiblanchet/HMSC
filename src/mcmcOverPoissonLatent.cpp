#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcOverPoissonLatent.h"

using namespace arma ;
using namespace Rcpp ;

// Gibbs sampling for a model that includes a set of latent variables only.
//' @rdname mcmcProbitX
//' @export
// [[Rcpp::export]]
RcppExport SEXP mcmcOverPoissonLatent(arma::mat& Y,
								 arma::mat& Ylatent,
								 arma::umat& Random,
								 arma::vec& residVar,
								 arma::field< arma::mat >& latent,
								 arma::field< arma::mat >& paramLatent,
								 arma::field< arma::mat >& shrinkLocal,
								 arma::field< arma::vec >& paramShrinkGlobal,
								 double priorResidVarScale,
								 double priorResidVarShape,
								 double priorShrinkLocal,
								 double priorShrinkOverallShape,
								 double priorShrinkOverallScale,
								 double priorShrinkSpeedShape,
								 double priorShrinkSpeedScale,
								 arma::vec& adapt,
								 arma::vec& redund,
								 int nRandom,
								 arma::vec& nRandomLev,
								 arma::vec& nLatent,
								 double nsp,
								 int nsite,
								 int niter,
								 int nburn,
								 int thin,
								 int verbose){

	// Define various objects
	mat EstModel = zeros<mat>(nsite,nsp);
	mat Yresid(nsite,nsp);

	// Define latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
	field<vec> shrinkGlobal(nRandom,1);
	field<mat> shrink(nRandom,1);

	// Initiate latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
	for (int i = 0; i < nRandom; i++) {
		shrinkGlobal(i,0) = cumprod(paramShrinkGlobal(i,0));

		mat shrinkLocalMat = trans(shrinkLocal(i,0));
		shrink(i,0) = trans(shrinkLocalMat.each_col() % shrinkGlobal(i,0));
	}

	// Define the result objects for burning
	field<mat> latentBurn(nburn/thin,nRandom);
	field<mat> paramLatentBurn(nburn/thin,nRandom);

	mat varDistBurn(nsp,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	field<mat> latentEst(nEst/thin,nRandom);
	field<mat> paramLatentEst(nEst/thin,nRandom);

	mat varDistEst(nsp,nEst/thin);

	// Define a burning counters
	int countBurn;
	countBurn = 0;

	// Define an estimation counters
	int countEst;
	countEst = 0;

	// Gibbs sampler
	for (int i = 0; i < niter; i++) {

		// Update paramLatent
		paramLatent = updateParamLatent(Ylatent, Random, residVar, paramLatent, latent, shrink, nRandom, nLatent, nsp, nsite);

		// Update latent
		latent = updateLatent(Ylatent, Random, residVar, paramLatent, latent, nRandom, nRandomLev, nLatent, nsp, nsite);

		//----------------
		// Sample Y latent
		//----------------
		// EstModel matrix full of zeros
		EstModel.zeros();

		// Add the effect of the latent variables
		for(int j = 0; j < nRandom; j++){
			mat latentMat = latent(j,0);
			EstModel = EstModel + latentMat.rows(Random.col(j))*trans(paramLatent(j,0));
		}

		Ylatent = sampleYlatentPoisson(Y, Ylatent, EstModel, residVar, nsp, nsite);

		//--------------------
		// Calculate residuals
		//--------------------
		Yresid = Ylatent-EstModel;

		//----------------
		// Update residVar
		//----------------
		residVar = updateResidVar(Yresid, residVar, priorResidVarScale, priorResidVarShape, nsp, nsite);

		//------------------------------
		// Shrinkage of latent variables
		//------------------------------
		for(int j = 0; j < nRandom; j++){
			// Update shrinkGlobal
			shrinkGlobal(j,0) = cumprod(paramShrinkGlobal(j,0));

			// Square each values in paramLatent(j,0)
			mat paramLatent2 = square(paramLatent(j,0));

			// Update shrinkLocal
			shrinkLocal(j,0) = updateShrinkLocal(shrinkLocal(j,0), priorShrinkLocal, shrinkGlobal(j,0), paramLatent2, nsp, nLatent(j));

			// Update paramShrinkGlobal
			paramShrinkGlobal(j,0) = updateParamShrinkGlobal(shrinkLocal(j,0), paramLatent2 ,paramShrinkGlobal(j,0), shrinkGlobal(j,0), priorShrinkOverallShape, priorShrinkOverallScale, priorShrinkSpeedShape, priorShrinkSpeedScale, nsp, nLatent(j));

			// Update shrink
			mat shrinkLocalMat = trans(shrinkLocal(j,0));
			shrink(j,0) = trans(shrinkLocalMat.each_col() % shrinkGlobal(j,0));
		}

		//-------------------------------------
		// Adapt the number of latent variables
		//-------------------------------------

		double probAdapt = 1/exp(adapt(0)+(adapt(1)*i));

		field<mat> tmpAdaptVar(7,1);

		for(int j = 0; j < nRandom; j++){
			tmpAdaptVar = adaptVar(paramLatent(j,0), latent(j,0), shrinkLocal(j,0), paramShrinkGlobal(j,0), shrinkGlobal(j,0), shrink(j,0), redund, priorShrinkLocal, priorShrinkSpeedShape, priorShrinkSpeedScale, probAdapt, nsp, nLatent(j), nRandomLev(j), i);

			latent(j,0) = tmpAdaptVar(0,0);
			nLatent(j) = tmpAdaptVar(1,0)(0,0);
			paramLatent(j,0) = tmpAdaptVar(2,0);
			shrinkLocal(j,0) = tmpAdaptVar(3,0);
			paramShrinkGlobal(j,0) = tmpAdaptVar(4,0);
			shrinkGlobal(j,0) = tmpAdaptVar(5,0);
			shrink(j,0) = tmpAdaptVar(6,0);
		}

		if(i<nburn && i%thin==0){
			// Save burning results
			for(int j = 0; j < nRandom; j++){
				paramLatentBurn(countBurn,j) = paramLatent(j,0);
				latentBurn(countBurn,j) = latent(j,0);
			}

			varDistBurn.col(countBurn) = residVar;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			for(int j = 0; j < nRandom; j++){
				paramLatentEst(countEst,j) = paramLatent(j,0);
				latentEst(countEst,j) = latent(j,0);
			}

			varDistEst.col(countEst) = residVar;

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
				Rcpp::List::create(Rcpp::Named("paramLatent", wrap(paramLatentBurn)),
								   Rcpp::Named("latent", wrap(latentBurn)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistBurn)))),
				Rcpp::List::create(Rcpp::Named("paramLatent", wrap(paramLatentEst)),
								   Rcpp::Named("latent", wrap(latentEst)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistEst)))));

}
