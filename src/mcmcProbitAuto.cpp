#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcProbitAuto.h"

using namespace arma ;
using namespace Rcpp ;

//[[Rcpp::export]]
RcppExport SEXP mcmcProbitAuto(arma::mat& Y,
							   arma::mat& Ylatent,
							   arma::field< arma::mat >& Auto,
							   arma::umat& RandomAuto,
							   arma::vec& residVar,
							   arma::field< arma::vec >& paramAuto,
							   arma::field< arma::mat >& latentAuto,
							   arma::field< arma::mat >& paramLatentAuto,
							   arma::field< arma::mat >& shrinkLocalAuto,
							   arma::field< arma::vec >& paramShrinkGlobalAuto,
							   double priorResidVarScale,
							   double priorResidVarShape,
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
							   int npriorParamAuto,
							   int niter,
							   int nburn,
							   int thin,
							   int verbose){

	// Define various objects
	mat EstModel = zeros<mat>(nsite,nsp);
	uvec Y0Loc = find(Y==0);
	uvec Y1Loc = find(Y==1);
	uvec YNALoc = find_nonfinite(Y);

	mat Yresid(nsite,nsp);
	double probAdapt;

	mat wAutoDet(npriorParamAuto,nAuto);
	field<cube> wAutoInv(nAuto,1);

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
	field<mat> latentAutoBurn(nburn/thin,nAuto);
	field<mat> paramLatentAutoBurn(nburn/thin,nAuto);
	field<vec> paramAutoBurn(nburn/thin,nAuto);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	field<mat> latentAutoEst(nEst/thin,nAuto);
	field<mat> paramLatentAutoEst(nEst/thin,nAuto);
	field<vec> paramAutoEst(nEst/thin,nAuto);

	// Define a burning counters
	int countBurn;
	countBurn = 0;

	// Define an estimation counters
	int countEst;
	countEst = 0;

	// Calculate Auto parameters that do not change in the Gibbs sampler
	field<cube> paramAutoNoGibb(2,1);
	for (int i = 0; i < nAuto; i++) {
		paramAutoNoGibb = fixParamAuto(Auto(i,0),priorParamAutoDist,nsp,nAutoLev,npriorParamAuto,i);
		wAutoDet.col(i) = vectorise(paramAutoNoGibb(0,0));
		wAutoInv(i,0) = paramAutoNoGibb(1,0);
	}

	// Convert NumericMatrix to arma::mat
	mat priorParamAutoDistArma = as<arma::mat>(priorParamAutoDist);

	// Gibbs sampler
	for (int i = 0; i < niter; i++) {
		//-----------------------
		// Update paramLatentAuto
		//-----------------------
		paramLatentAuto = updateParamLatent(Ylatent, RandomAuto, residVar, paramLatentAuto, latentAuto, shrinkAuto, nAuto, nLatentAuto, nsp, nsite);

		//------------------
		// Update latentAuto
		//------------------
		latentAuto = updateLatentAuto(Ylatent, RandomAuto, residVar, paramAuto, wAutoInv, paramLatentAuto, latentAuto, priorParamAutoDistArma, nAuto, nAutoLev, nLatentAuto, nsp, nsite);

		//-----------------
		// Update paramAuto
		//-----------------
		paramAuto = updateParamAuto(latentAuto, paramAuto,npriorParamAuto,wAutoDet,wAutoInv, priorParamAutoWeight, priorParamAutoDist,nLatentAuto,nAuto);

		//----------------
		// Sample Y latent
		//----------------
		// EstModel matrix full of zeros
		EstModel.zeros();
		// Add the effect of the autocorrelated latent variables
		for(int j = 0; j < nAuto; j++){
			mat latentAutoMat = latentAuto(j,0);
			EstModel = EstModel + (latentAutoMat.rows(RandomAuto.col(j))*trans(paramLatentAuto(j,0)));
		}

		Ylatent = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, Ylatent, EstModel, residVar, nsp, nsite);
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

		//-------------
		// Save results
		//-------------
		if(i<nburn && i%thin==0){
			// Save burning results
			for(int j = 0; j < nAuto; j++){
				paramLatentAutoBurn(countBurn,j) = paramLatentAuto(j,0);
				latentAutoBurn(countBurn,j) = latentAuto(j,0);
				paramAutoBurn(countBurn,j) = paramAuto(j,0);
			}

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			for(int j = 0; j < nAuto; j++){
				paramLatentAutoEst(countEst,j) = paramLatentAuto(j,0);
				latentAutoEst(countEst,j) = latentAuto(j,0);
				paramAutoEst(countEst,j) = paramAuto(j,0);
			}

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
				Rcpp::List::create(Rcpp::Named("paramLatentAuto", wrap(paramLatentAutoBurn)),
								   Rcpp::Named("latentAuto", wrap(latentAutoBurn)),
								   Rcpp::Named("paramAuto", wrap(paramAutoBurn))),
				Rcpp::List::create(Rcpp::Named("paramLatentAuto", wrap(paramLatentAutoEst)),
								   Rcpp::Named("latentAuto", wrap(latentAutoEst)),
								   Rcpp::Named("paramAuto", wrap(paramAutoEst))));

}
