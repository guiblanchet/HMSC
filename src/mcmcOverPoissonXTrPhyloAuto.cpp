#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcOverPoissonXTrPhyloAuto.h"

using namespace arma ;
using namespace Rcpp ;

//' @rdname mcmcProbitX
//' @export
//[[Rcpp::export]]
RcppExport SEXP mcmcOverPoissonXTrPhyloAuto(arma::mat& Y,
									   arma::mat& Ylatent,
									   arma::mat& X,
									   arma::mat& Tr,
									   arma::mat& Phylo,
									   arma::mat& iPhylo,
									   arma::field< arma::mat >& Auto,
									   arma::umat& RandomAuto,
									   arma::mat& paramX,
									   arma::mat& paramTr,
									   arma::mat& precX,
									   double paramPhylo,
									   arma::vec& residVar,
									   arma::field< arma::vec >& paramAuto,
									   arma::field< arma::mat >& latentAuto,
									   arma::field< arma::mat >& paramLatentAuto,
									   arma::field< arma::mat >& shrinkLocalAuto,
									   arma::field< arma::vec >& paramShrinkGlobalAuto,
									   arma::mat& priorParamTr,
									   arma::mat& priorVarTr,
									   arma::mat& priorVarXScaleMat,
									   double priorVarXDf,
									   double priorResidVarScale,
									   double priorResidVarShape,
									   arma::mat& priorParamPhyloWeight,
									   Rcpp::NumericVector& priorParamPhyloGrid,
									   arma::mat& priorParamAutoWeight,
									   Rcpp::NumericMatrix& priorParamAutoDist,
									   double priorShrinkLocal,
									   double priorShrinkOverallShape,
									   double priorShrinkOverallScale,
									   double priorShrinkSpeedShape,
									   double priorShrinkSpeedScale,
									   arma::vec& adapt,
									   arma::vec& redund,
									   int nAuto,
									   arma::vec& nAutoLev,
									   arma::vec& nLatentAuto,
									   double nsp,
									   int nsite,
									   int nparamX,
									   int nTr,
									   int nparamPhylo,
									   int npriorParamAuto,
									   int niter,
									   int nburn,
									   int thin,
									   int verbose){

	// Define various objects
	mat EstModel = zeros<mat>(nsite,nsp);
	mat Yresid(nsite,nsp);
	mat varX(nparamX,nparamX);
	double probAdapt;

	mat wAutoDet(npriorParamAuto,nAuto);
	field<cube> wAutoInv(nAuto,1);

	int paramPhyloPointer;
	mat wPhyloInvMat(nsp,nsp);
	vec wPhyloDet(nparamPhylo);
	cube wPhyloInv(nsp,nsp,nparamPhylo);

	// Define meansParamX
	mat meansParamX(nsp,nparamX);
	meansParamX = trans(paramTr*Tr);

	///////////////////////////////
	// Autocorrelated random effect
	///////////////////////////////
	// Define latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
	field<vec> shrinkGlobalAuto(nAuto,1);
	field<mat> shrinkAuto(nAuto,1);

	// Initiate latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
	for (int i = 0; i < nAuto; i++) {
		shrinkGlobalAuto(i,0) = cumprod(paramShrinkGlobalAuto(i,0));

		mat shrinkLocalAutoMat = trans(shrinkLocalAuto(i,0));
		shrinkAuto(i,0) = trans(shrinkLocalAutoMat.each_col() % shrinkGlobalAuto(i,0));
	}

	/////////////////
	// Result objects
	/////////////////

	// Define the result objects for burning
	cube paramTrBurn(nparamX,nTr,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);

	vec paramPhyloBurn(nburn/thin);

	field<mat> latentAutoBurn(nburn/thin,nAuto);
	field<mat> paramLatentAutoBurn(nburn/thin,nAuto);
	field<vec> paramAutoBurn(nburn/thin,nAuto);

	mat varDistBurn(nsp,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	cube paramTrEst(nparamX,nTr,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);

	vec paramPhyloEst(nEst/thin);

	field<mat> latentAutoEst(nEst/thin,nAuto);
	field<mat> paramLatentAutoEst(nEst/thin,nAuto);
	field<vec> paramAutoEst(nEst/thin,nAuto);

	mat varDistEst(nsp,nEst/thin);

	// Calculate Phylo parameters that do not change in the Gibbs sampler
	field<cube> paramPhyloNoGibb = fixParamPhylo(Phylo,iPhylo,paramPhylo,priorParamPhyloGrid,nsp,nparamPhylo);
	wPhyloDet = vectorise(paramPhyloNoGibb(0,0));
	wPhyloInv = paramPhyloNoGibb(1,0);
	mat priorParamPhyloGridArma = as<arma::vec>(priorParamPhyloGrid); // Convert NumericVector to arma::vec

	// Calculate Auto parameters that do not change in the Gibbs sampler
	field<cube> paramAutoNoGibb(2,1);
	for (int i = 0; i < nAuto; i++) {
		paramAutoNoGibb = fixParamAuto(Auto(i,0),priorParamAutoDist,nsp,nAutoLev,npriorParamAuto,i);
		wAutoDet.col(i) = vectorise(paramAutoNoGibb(0,0));
		wAutoInv(i,0) = paramAutoNoGibb(1,0);
	}

	// Convert NumericMatrix to arma::mat
	mat priorParamAutoDistArma = as<arma::mat>(priorParamAutoDist);

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

		// Remove influence of X variables
		Yresid = Ylatent-EstModel;

		//-----------------------
		// Update paramLatentAuto
		//-----------------------
		paramLatentAuto = updateParamLatent(Yresid, RandomAuto, residVar, paramLatentAuto, latentAuto, shrinkAuto, nAuto, nLatentAuto, nsp, nsite);

		//------------------
		// Update latentAuto
		//------------------
		latentAuto = updateLatentAuto(Yresid, RandomAuto, residVar, paramAuto, wAutoInv, paramLatentAuto, latentAuto, priorParamAutoDistArma, nAuto, nAutoLev, nLatentAuto, nsp, nsite);

		//-----------------
		// Update paramAuto
		//-----------------
		paramAuto = updateParamAuto(latentAuto, paramAuto,npriorParamAuto,wAutoDet,wAutoInv, priorParamAutoWeight, priorParamAutoDist,nLatentAuto,nAuto);

		//----------------
		// Sample Y latent
		//----------------
		// Add the effect of the autocorrelated latent variables
		for(int j = 0; j < nAuto; j++){
			mat latentAutoMat = latentAuto(j,0);
			EstModel = EstModel + (latentAutoMat.rows(RandomAuto.col(j))*diagmat(paramAuto(j,0))*trans(paramLatentAuto(j,0))); // Not sure if multiplying by paramAuto is the way to go.
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

		//--------------
		// Update paramX
		//--------------
		// Remove the effect of the latent variables
		Yresid = Ylatent;
		for(int j = 0; j < nAuto; j++){
			mat latentAutoMat = latentAuto(j,0);
			Yresid = Yresid - (latentAutoMat.rows(RandomAuto.col(j))*diagmat(paramAuto(j,0))*trans(paramLatentAuto(j,0)));
		}

		meansParamX = trans(paramTr*Tr);
		paramX = updateParamXPhylo(Yresid, X, paramX, meansParamX, precX, residVar, wPhyloInvMat, nsp, nsite, nparamX);

		//-------------
		// Update precX
		//-------------
		precX = updatePrecXPhylo(meansParamX, paramX, precX, wPhyloInvMat, priorVarXScaleMat, priorVarXDf, nsp);
		varX = precX.i();

		//---------------
		// Update paramTr
		//---------------
		paramTr = updateParamTrPhylo(Tr, paramX, paramTr, precX, wPhyloInvMat, priorParamTr, priorVarTr, nTr, nparamX);

		//------------------
		// Update paramPhylo
		//------------------
		paramPhylo = updateParamPhylo(paramX,meansParamX,paramPhylo,precX,wPhyloDet,wPhyloInv,nparamX,nparamPhylo,nsp,priorParamPhyloWeight,priorParamPhyloGrid);

		//---------------------------------------------
		// Shrinkage of autocorrelated latent variables
		//---------------------------------------------
		for(int j = 0; j < nAuto; j++){
			// Update shrinkLocal
			mat paramLatentAuto2 = square(paramLatentAuto(j,0));
			shrinkLocalAuto(j,0) = updateShrinkLocal(shrinkLocalAuto(j,0), priorShrinkLocal, shrinkGlobalAuto(j,0), paramLatentAuto2, nsp, nLatentAuto(j));

			// Update paramShrinkGlobal
			paramShrinkGlobalAuto(j,0) = updateParamShrinkGlobal(shrinkLocalAuto(j,0), paramLatentAuto2 ,paramShrinkGlobalAuto(j,0), shrinkGlobalAuto(j,0), priorShrinkOverallShape, priorShrinkOverallScale, priorShrinkSpeedShape, priorShrinkSpeedScale, nsp, nLatentAuto(j));

			// Update shrinkGlobal
			shrinkGlobalAuto(j,0) = cumprod(paramShrinkGlobalAuto(j,0));

			// Update shrink
			mat shrinkLocalAutoMat = trans(shrinkLocalAuto(j,0));
			shrinkAuto(j,0) = trans(shrinkLocalAutoMat.each_col() % shrinkGlobalAuto(j,0));
		}

		//----------------------------------------------------
		// Adapt the number of autocorrelated latent variables
		//----------------------------------------------------
		probAdapt = 1/exp(adapt(0)+(adapt(1)*i));

		field<mat> tmpAdaptVarAuto(8,1);

		for(int j = 0; j < nAuto; j++){
			tmpAdaptVarAuto = adaptVarAuto(paramLatentAuto(j,0), latentAuto(j,0),paramAuto(j,0), shrinkLocalAuto(j,0), paramShrinkGlobalAuto(j,0),  shrinkGlobalAuto(j,0), shrinkAuto(j,0), redund, priorShrinkLocal, priorShrinkSpeedShape, priorShrinkSpeedScale, priorParamAutoDistArma, probAdapt, nsp, nLatentAuto(j), nAutoLev(j), i);

			latentAuto(j,0) = tmpAdaptVarAuto(0,0);
			nLatentAuto(j) = tmpAdaptVarAuto(1,0)(0,0);
			paramLatentAuto(j,0) = tmpAdaptVarAuto(2,0);
			shrinkLocalAuto(j,0) = tmpAdaptVarAuto(3,0);
			paramShrinkGlobalAuto(j,0) = tmpAdaptVarAuto(4,0);
			shrinkGlobalAuto(j,0) = tmpAdaptVarAuto(5,0);
			shrinkAuto(j,0) = tmpAdaptVarAuto(6,0);
			paramAuto(j,0) = tmpAdaptVarAuto(7,0);
		}

		if(i<nburn && i%thin==0){
			// Save burning results
			paramTrBurn.slice(countBurn) = paramTr;
			varXBurn.slice(countBurn) = precX.i();
			paramXBurn.slice(countBurn) = paramX;
			paramPhyloBurn(countBurn) = paramPhylo;

			for(int j = 0; j < nAuto; j++){
				paramLatentAutoBurn(countBurn,j) = paramLatentAuto(j,0);
				latentAutoBurn(countBurn,j) = latentAuto(j,0);
				paramAutoBurn(countBurn,j) = paramAuto(j,0);
			}

			varDistBurn.col(countBurn) = residVar;

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			paramTrEst.slice(countEst) = paramTr;
			varXEst.slice(countEst) = precX.i();
			paramXEst.slice(countEst) = paramX;
			paramPhyloEst(countEst) = paramPhylo;

			for(int j = 0; j < nAuto; j++){
				paramLatentAutoEst(countEst,j) = paramLatentAuto(j,0);
				latentAutoEst(countEst,j) = latentAuto(j,0);
				paramAutoEst(countEst,j) = paramAuto(j,0);
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
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXBurn)),
								   Rcpp::Named("paramTr", wrap(paramTrBurn)),
								   Rcpp::Named("varX", wrap(varXBurn)),
								   Rcpp::Named("paramPhylo", wrap(paramPhyloBurn)),
								   Rcpp::Named("paramLatentAuto", wrap(paramLatentAutoBurn)),
								   Rcpp::Named("latentAuto", wrap(latentAutoBurn)),
								   Rcpp::Named("paramAuto", wrap(paramAutoBurn)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistBurn)))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("paramTr", wrap(paramTrEst)),
								   Rcpp::Named("varX", wrap(varXEst)),
								   Rcpp::Named("paramPhylo", wrap(paramPhyloEst)),
								   Rcpp::Named("paramLatentAuto", wrap(paramLatentAutoEst)),
								   Rcpp::Named("latentAuto", wrap(latentAutoEst)),
								   Rcpp::Named("paramAuto", wrap(paramAutoEst)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistEst)))));

}
