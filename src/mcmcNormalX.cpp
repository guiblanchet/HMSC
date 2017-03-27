#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcNormalX.h"

using namespace arma ;
using namespace Rcpp ;

//' @rdname mcmcProbitX
//' @export
//[[Rcpp::export]]
RcppExport SEXP mcmcNormalX(arma::mat& Ylatent,
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
	mat varDistBurn(nsp,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	mat meansParamXEst(nparamX,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);
	mat varDistEst(nsp,nEst/thin);

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

		// Calculate residuals
		Yresid = Ylatent-EstModel;

		// Update residVar
		residVar = updateResidVar(Yresid, residVar, priorResidVarScale, priorResidVarShape, nsp, nsite);

		// Update paramX
		meansParamXRep = trans(repmat(meansParamX,1,nsp));
		paramX = updateParamX1(Ylatent,X,meansParamXRep,precX, paramX, residVar, nsp, nsite, nparamX);

		// Update precX
		precX = updatePrecX(meansParamX,priorVarXScaleMat, priorVarXDf, paramX, precX, nsp);

		// Update meanparamX
		meansParamX = updateMeansParamX(priorMeansParamX, priorVarXScaleMat, priorVarXDf, paramX, meansParamX, precX, nsp, nparamX);

		if(i<nburn && i%thin==0){
			// Save burning results
			meansParamXBurn.col(countBurn) = meansParamX;
			varDistBurn.col(countBurn) = residVar;
			varXBurn.slice(countBurn) = precX.i();
			paramXBurn.slice(countBurn) = paramX;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			meansParamXEst.col(countEst) = meansParamX;
			varDistEst.col(countEst) = residVar;
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
								   Rcpp::Named("varX", wrap(varXBurn)),
								   Rcpp::Named("varNormal", wrap(trans(1/varDistBurn)))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("meansParamX", wrap(trans(meansParamXEst))),
								   Rcpp::Named("varX", wrap(varXEst)),
								   Rcpp::Named("varNormal", wrap(trans(1/varDistEst)))));

}
