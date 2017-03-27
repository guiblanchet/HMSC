#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcPoissonXPhylo.h"

using namespace arma ;
using namespace Rcpp ;

//' @rdname mcmcProbitX
//' @export
//[[Rcpp::export]]
RcppExport SEXP mcmcPoissonXPhylo(arma::mat& Y,
								 arma::mat& Ylatent,
								 arma::mat& X,
								 arma::mat& Phylo,
								 arma::mat& iPhylo,
								 arma::mat& paramX,
								 arma::mat& meansParamX,
								 double paramPhylo,
								 arma::mat& precX,
								 arma::vec& residVar,
								 arma::mat& priorMeansParamX,
								 arma::mat& priorVarMeansParamX,
								 arma::mat& priorVarXScaleMat,
								 double priorVarXDf,
								 double priorResidVarScale,
								 double priorResidVarShape,
								 arma::mat& priorParamPhyloWeight,
								 Rcpp::NumericVector& priorParamPhyloGrid,
								 double nsp,
								 int nsite,
								 int nparamX,
								 int nparamPhylo,
								 int niter,
								 int nburn,
								 int thin,
								 int verbose){

	// Define various objects
	mat EstModel(nsite,nsp);
	mat eyeSp(nsp,nsp,fill::eye);
	mat varX(nparamX,nparamX);
	int paramPhyloPointer;
	mat wPhyloInvMat(nsp,nsp);

	vec wPhyloDet(nparamPhylo);
	cube wPhyloInv(nsp,nsp,nparamPhylo);

	// Define the result objects for burning
	mat meansParamXBurn(nparamX,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);
	vec paramPhyloBurn(nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	mat meansParamXEst(nparamX,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);
	vec paramPhyloEst(nEst/thin);

	// Calculate Phylo parameters that do not change in the Gibbs sampler
	field<cube> paramPhyloNoGibb = fixParamPhylo(Phylo,iPhylo,paramPhylo,priorParamPhyloGrid,nsp,nparamPhylo);
	wPhyloDet = vectorise(paramPhyloNoGibb(0,0));
	wPhyloInv = paramPhyloNoGibb(1,0);
	mat priorParamPhyloGridArma = as<arma::vec>(priorParamPhyloGrid); // Convert NumericVector to arma::vec

	// Define a burning counters
	int countBurn;
	countBurn = 0;

	// Define an estimation counters
	int countEst;
	countEst = 0;

	// Gibbs sampler
	for (int i = 0; i < niter; i++) {
		// Extract the matrix in wPhyloInv used to estimate paramX and paramTr
		paramPhyloPointer = as_scalar(find(priorParamPhyloGridArma==paramPhylo)); // Find which wPhyloInv in wPhyloInv is the one used
		wPhyloInvMat = wPhyloInv.slice(paramPhyloPointer); // Extract the right part of wPhyloInv

		// Calculate the model estimate
		EstModel = X*trans(paramX);

		// Sample Y latent
		Ylatent = sampleYlatentPoisson(Y, Ylatent, EstModel, residVar, nsp, nsite);

		// Update paramX
		mat meansParamXRep = trans(repmat(meansParamX,1,nsp));
		paramX = updateParamXPhylo(Ylatent, X, paramX, meansParamXRep, precX, residVar, wPhyloInvMat, nsp, nsite, nparamX);

		// Update precX
		precX = updatePrecXPhylo(meansParamXRep, paramX, precX, wPhyloInvMat, priorVarXScaleMat, priorVarXDf, nsp);
		varX = precX.i();

		// Update meanparamX
		meansParamX = updateMeansParamXPhylo(paramX, meansParamX, precX, wPhyloInvMat, priorMeansParamX, priorVarMeansParamX, nsp);

		// Update paramPhylo
		paramPhylo = updateParamPhylo(paramX,meansParamXRep,paramPhylo,precX,wPhyloDet,wPhyloInv,nparamX,nparamPhylo,nsp,priorParamPhyloWeight,priorParamPhyloGrid);

		if(i<nburn && i%thin==0){
			// Save burning results
			meansParamXBurn.col(countBurn) = meansParamX;
			varXBurn.slice(countBurn) = varX;
			paramXBurn.slice(countBurn) = paramX;
			paramPhyloBurn(countBurn) = paramPhylo;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			meansParamXEst.col(countEst) = meansParamX;
			varXEst.slice(countEst) = varX;
			paramXEst.slice(countEst) = paramX;
			paramPhyloEst(countEst) = paramPhylo;

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
								   Rcpp::Named("paramPhylo", wrap(paramPhyloBurn))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("meansParamX", wrap(trans(meansParamXEst))),
								   Rcpp::Named("varX", wrap(varXEst)),
								   Rcpp::Named("paramPhylo", wrap(paramPhyloEst))));

}
