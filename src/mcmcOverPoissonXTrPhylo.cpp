#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcOverPoissonXTrPhylo.h"

using namespace arma ;
using namespace Rcpp ;

//' @rdname mcmcProbitX
//' @export
//[[Rcpp::export]]
RcppExport SEXP mcmcOverPoissonXTrPhylo(arma::mat& Y,
								   arma::mat& Ylatent,
								   arma::mat& X,
								   arma::mat& Tr,
								   arma::mat& Phylo,
								   arma::mat& iPhylo,
								   arma::mat& paramX,
								   arma::mat& paramTr,
								   double paramPhylo,
								   arma::mat& precX,
								   arma::vec& residVar,
								   arma::mat& priorParamTr,
								   arma::mat& priorVarTr,
								   arma::mat& priorVarXScaleMat,
								   double priorVarXDf,
								   double priorResidVarScale,
								   double priorResidVarShape,
								   arma::mat& priorParamPhyloWeight,
								   Rcpp::NumericVector& priorParamPhyloGrid,
								   double nsp,
								   int nsite,
								   int nparamX,
								   int nTr,
								   int nparamPhylo,
								   int niter,
								   int nburn,
								   int thin,
								   int verbose){

	// Define various objects
	mat EstModel(nsite,nsp);
	mat Yresid(nsite,nsp);
	mat varX(nparamX,nparamX);
	mat eyeSp(nsp,nsp,fill::eye);
	int paramPhyloPointer;
	mat wPhyloInvMat(nsp,nsp);

	vec wPhyloDet(nparamPhylo);
	cube wPhyloInv(nsp,nsp,nparamPhylo);

	// Define meansParamX
	mat meansParamX(nsp,nparamX);
	meansParamX = trans(paramTr*Tr);

	// Define the result objects for burning
	cube paramTrBurn(nparamX,nTr,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);
	vec paramPhyloBurn(nburn/thin);

	mat varDistBurn(nsp,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	cube paramTrEst(nparamX,nTr,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);
	vec paramPhyloEst(nEst/thin);

	mat varDistEst(nsp,nEst/thin);

	// Calculate Phylo parameters that do not change in the Gibbs sampler
	field<cube> ParamPhyloNoGibb = fixParamPhylo(Phylo,iPhylo,paramPhylo,priorParamPhyloGrid,nsp,nparamPhylo);
	wPhyloDet = vectorise(ParamPhyloNoGibb(0,0));
	wPhyloInv = ParamPhyloNoGibb(1,0);

	// Define a burning counters
	int countBurn;
	countBurn = 0;

	// Define an estimation counters
	int countEst;
	countEst = 0;

	// Gibbs sampler
	for (int i = 0; i < niter; i++) {
		// Extract the matrix in wPhyloInv used to estimate paramX and paramTr
		mat priorParamPhyloGridArma = as<arma::vec>(priorParamPhyloGrid); // Convert NumericVector to arma::vec
		paramPhyloPointer = as_scalar(find(priorParamPhyloGridArma==paramPhylo)); // Find which wPhyloInv in wPhyloInv is the one used
		wPhyloInvMat = wPhyloInv.slice(paramPhyloPointer); // Extract the right part of wPhyloInv

		// Recalculate meansParamX
		meansParamX = trans(paramTr*Tr);

		// Calculate the model estimate
		EstModel = X*trans(paramX);

		// Sample Y latent
		Ylatent = sampleYlatentPoisson(Y, Ylatent, EstModel, residVar, nsp, nsite);

		// Calculate residuals
		Yresid = Ylatent-EstModel;

		// Update residVar
		residVar = updateResidVar(Yresid, residVar, priorResidVarScale, priorResidVarShape, nsp, nsite);

		// Update paramX
		paramX = updateParamXPhylo(Ylatent, X, paramX, meansParamX, precX, residVar, wPhyloInvMat, nsp, nsite, nparamX);

		// Update precX and calculate varX
		precX = updatePrecXPhylo(meansParamX, paramX, precX, wPhyloInvMat, priorVarXScaleMat, priorVarXDf, nsp);
		varX = precX.i();

		// Update paramTr
		paramTr = updateParamTrPhylo(Tr, paramX, paramTr, precX, wPhyloInvMat, priorParamTr, priorVarTr, nTr, nparamX);

		// Update paramPhylo
		paramPhylo = updateParamPhylo(paramX,meansParamX,paramPhylo,precX,wPhyloDet,wPhyloInv,nparamX,nparamPhylo,nsp,priorParamPhyloWeight,priorParamPhyloGrid);

		if(i<nburn && i%thin==0){
			// Save burning results
			paramTrBurn.slice(countBurn) = paramTr;
			varXBurn.slice(countBurn) = varX;
			paramXBurn.slice(countBurn) = paramX;
			paramPhyloBurn(countBurn) = paramPhylo;

			varDistBurn.col(countBurn) = residVar;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			paramTrEst.slice(countEst) = paramTr;
			varXEst.slice(countEst) = varX;
			paramXEst.slice(countEst) = paramX;
			paramPhyloEst(countEst) = paramPhylo;

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
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXBurn)),
								   Rcpp::Named("paramTr", wrap(paramTrBurn)),
								   Rcpp::Named("varX", wrap(varXBurn)),
								   Rcpp::Named("paramPhylo", wrap(paramPhyloBurn)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistBurn)))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("paramTr", wrap(paramTrEst)),
								   Rcpp::Named("varX", wrap(varXEst)),
								   Rcpp::Named("paramPhylo", wrap(paramPhyloEst)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistEst)))));

}
