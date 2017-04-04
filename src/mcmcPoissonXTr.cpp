#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcPoissonXTr.h"

using namespace arma ;
using namespace Rcpp ;

// Gibbs sampling for a model that includes a set of explanatory variable only.
//[[Rcpp::export]]
RcppExport SEXP mcmcPoissonXTr(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& X,
							  arma::mat& Tr,
							  arma::mat& paramX,
							  arma::mat& paramTr,
							  arma::mat& precX,
							  arma::vec& residVar,
							  arma::mat& priorParamTr,
							  arma::mat& priorVarTr,
							  arma::mat& priorVarXScaleMat,
							  double priorVarXDf,
							  double priorResidVarScale,
							  double priorResidVarShape,
							  double nsp,
							  int nsite,
							  int nparamX,
							  int nTr,
							  int niter,
							  int nburn,
							  int thin,
							  int verbose){

	// Define various objects
	mat EstModel(nsite,nsp);
	mat varX(nparamX,nparamX);

	// Define meansParamX
	mat meansParamX(nsp,nparamX);
	meansParamX = trans(paramTr*Tr);

	// Define the result objects for burning
	cube paramTrBurn(nparamX,nTr,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	cube paramTrEst(nparamX,nTr,nEst/thin);
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

		// Sample Y latent
		Ylatent = sampleYlatentPoisson(Y, Ylatent, EstModel, residVar, nsp, nsite);

		// Update paramX
		paramX = updateParamX(Ylatent,X,meansParamX,precX, paramX, nsp, nsite, nparamX);

		// Update precX and calculate varX
		precX = updatePrecXTr(Tr, priorVarXScaleMat, priorVarXDf, paramTr, paramX, precX, nsp, nparamX);
		varX = precX.i();

		// Update paramTr
		paramTr = updateParamTr(Tr, paramX, paramTr, precX, priorVarTr, priorParamTr, nparamX, nTr);

		// Recalculate meansParamX
		meansParamX = trans(paramTr*Tr);

		if(i<nburn && i%thin==0){
			// Save burning results
			paramTrBurn.slice(countBurn) = paramTr;
			varXBurn.slice(countBurn) = varX;
			paramXBurn.slice(countBurn) = paramX;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			paramTrEst.slice(countEst) = paramTr;
			varXEst.slice(countEst) = varX;
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
								   Rcpp::Named("paramTr", wrap(paramTrBurn)),
								   Rcpp::Named("varX", wrap(varXBurn))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("paramTr", wrap(paramTrEst)),
								   Rcpp::Named("varX", wrap(varXEst))));

}
