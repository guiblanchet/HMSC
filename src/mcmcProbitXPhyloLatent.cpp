#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcProbitXPhyloLatent.h"

using namespace arma ;
using namespace Rcpp ;

// Gibbs sampling for a model that includes a set of explanatory variable, latent variables and phylogeny.
// [[Rcpp::export]]
RcppExport SEXP mcmcProbitXPhyloLatent(arma::mat& Y,
									   arma::mat& Ylatent,
									   arma::mat& X,
									   arma::mat& Phylo,
									   arma::mat& iPhylo,
									   arma::umat& Random,
									   arma::mat& paramX,
									   arma::mat& meansParamX,
									   arma::mat& precX,
									   double paramPhylo,
									   arma::vec& residVar,
									   arma::mat& priorMeansParamX,
									   arma::mat& priorVarMeansParamX,
									   arma::mat& priorVarXScaleMat,
									   double priorVarXDf,
									   double priorResidVarScale,
									   double priorResidVarShape,
									   arma::mat& priorParamPhyloWeight,
									   Rcpp::NumericVector& priorParamPhyloGrid,
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
									   int nparamX,
									   int nparamPhylo,
									   int niter,
									   int nburn,
									   int thin,
									   int verbose){

	// Define various objects
	mat EstModel(nsite,nsp);
	uvec Y0Loc = find(Y==0);
	uvec Y1Loc = find(Y==1);
	uvec YNALoc = find_nonfinite(Y);

	mat Yresid(nsite,nsp);
	mat varX(nparamX,nparamX);
	int paramPhyloPointer;
	mat wPhyloInvMat(nsp,nsp);

	vec wPhyloDet(nparamPhylo);
	cube wPhyloInv(nsp,nsp,nparamPhylo);

	// Define latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
	field<mat> latent(nRandom,1);
	field<mat> paramLatent(nRandom,1);
	field<mat> shrinkLocal(nRandom,1);
	field<vec> paramShrinkGlobal(nRandom,1);
	field<vec> shrinkGlobal(nRandom,1);
	field<mat> shrink(nRandom,1);

	// Initiate latent variables, their parameters, and the different parameters associated to the shrinkage of the latent variables
	for (int i = 0; i < nRandom; i++) {
		latent(i,0) = randn(nRandomLev(i),nLatent(i));
		paramLatent(i,0) = randn(nsp,nLatent(i));
		shrinkLocal(i,0) = randg(nsp,nLatent(i),distr_param(priorShrinkLocal/2,2/priorShrinkLocal));
		paramShrinkGlobal(i,0) = join_cols(randg(1,distr_param(priorShrinkOverallShape,1/priorShrinkOverallScale)),randg(nLatent(i)-1,distr_param(priorShrinkSpeedShape,1/priorShrinkSpeedScale)));
		shrinkGlobal(i,0) = cumprod(paramShrinkGlobal(i,0));

		mat shrinkLocalMat = trans(shrinkLocal(i,0));
		shrink(i,0) = trans(shrinkLocalMat.each_col() % shrinkGlobal(i,0));
	}

	// Define the result objects for burning
	mat meansParamXBurn(nparamX,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);
	vec paramPhyloBurn(nburn/thin);
	field<mat> latentBurn(nburn/thin,nRandom);
	field<mat> paramLatentBurn(nburn/thin,nRandom);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	mat meansParamXEst(nparamX,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);
	vec paramPhyloEst(nEst/thin);
	field<mat> latentEst(nEst/thin,nRandom);
	field<mat> paramLatentEst(nEst/thin,nRandom);

	// Define a burning counters
	int countBurn;
	countBurn = 0;

	// Define an estimation counters
	int countEst;
	countEst = 0;

	// Calculate Phylo parameters that do not change in the Gibbs sampler
	// Inverse phylogeny
	field<cube> ParamPhyloNoGibb = fixParamPhylo(Phylo,iPhylo,paramPhylo,priorParamPhyloGrid,nsp,nparamPhylo);
	wPhyloDet = vectorise(ParamPhyloNoGibb(0,0));
	wPhyloInv = ParamPhyloNoGibb(1,0);

	// Gibbs sampler
	for (int i = 0; i < niter; i++) {
		// Extract the matrix in wPhyloInv used to estimate paramX and paramTr
		mat priorParamPhyloGridArma = as<arma::vec>(priorParamPhyloGrid); // Convert NumericVector to arma::vec
		paramPhyloPointer = as_scalar(find(priorParamPhyloGridArma==paramPhylo)); // Find which wPhyloInv in wPhyloInv is the one used
		wPhyloInvMat = wPhyloInv.slice(paramPhyloPointer); // Extract the right part of wPhyloInv

		// Calculate the model estimate
		EstModel = X*trans(paramX);

		// Remove influence of X variables
		Yresid = Ylatent-EstModel;

		// Update paramLatent
		paramLatent = updateParamLatent(Yresid, Random, residVar, paramLatent, latent, shrink, nRandom, nLatent, nsp, nsite);

		// Update latent
		latent = updateLatent(Yresid, Random, residVar, paramLatent, latent, nRandom, nRandomLev, nLatent, nsp, nsite);

		//----------------
		// Sample Y latent
		//----------------
		// Remove the effect of the latent variables
		for(int j = 0; j < nRandom; j++){
			mat latentMat = latent(j,0);
			EstModel = EstModel + latentMat.rows(Random.col(j))*trans(paramLatent(j,0));
		}

		Ylatent = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, Ylatent, EstModel, residVar, nsp, nsite);

		//--------------
		// Update paramX
		//--------------
		// Remove the effect of the latent variables
		Yresid = Ylatent;
		for(int j = 0; j < nRandom; j++){
			mat latentMat = latent(j,0);
			Yresid = Yresid - latentMat.rows(Random.col(j))*trans(paramLatent(j,0));
		}

		mat meansParamXRep = trans(repmat(meansParamX,1,nsp));
		paramX = updateParamXPhylo(Yresid, X, paramX, meansParamXRep, precX, residVar, wPhyloInvMat, nsp, nsite, nparamX);

		//-------------
		// Update precX
		//-------------
		precX = updatePrecXPhylo(meansParamXRep, paramX, precX, wPhyloInvMat, priorVarXScaleMat, priorVarXDf, nsp);
		varX = precX.i();

		//------------------
		// Update meanparamX
		//------------------
		meansParamX = updateMeansParamXPhylo(paramX, meansParamX, precX, wPhyloInvMat, priorMeansParamX, priorVarMeansParamX, nsp);

		//------------------
		// Update paramPhylo
		//------------------
		paramPhylo = updateParamPhylo(paramX,meansParamXRep,paramPhylo,precX,wPhyloDet,wPhyloInv,nparamX,nparamPhylo,nsp,priorParamPhyloWeight,priorParamPhyloGrid);

		//------------------------------
		// Shrinkage of latent variables
		//------------------------------
		for(int j = 0; j < nRandom; j++){
			// Update shrinkLocal

			mat paramLatent2 = square(paramLatent(j,0));
			shrinkLocal(j,0) = updateShrinkLocal(shrinkLocal(j,0), priorShrinkLocal, shrinkGlobal(j,0), paramLatent2, nsp, nLatent(j));

			// Update paramShrinkGlobal
			paramShrinkGlobal(j,0) = updateParamShrinkGlobal(shrinkLocal(j,0), paramLatent2 ,paramShrinkGlobal(j,0), shrinkGlobal(j,0), priorShrinkOverallShape, priorShrinkOverallScale, priorShrinkSpeedShape, priorShrinkSpeedScale, nsp, nLatent(j));

			// Update shrinkGlobal
			shrinkGlobal(j,0) = cumprod(paramShrinkGlobal(j,0));

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
			tmpAdaptVar = adaptVar(paramLatent(j,0), latent(j,0), shrinkLocal(j,0), paramShrinkGlobal(j,0),  shrinkGlobal(j,0), shrink(j,0), redund, priorShrinkLocal, priorShrinkSpeedShape, priorShrinkSpeedScale, probAdapt, nsp, nLatent(j), nRandomLev(j), i);

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
			meansParamXBurn.col(countBurn) = meansParamX;
			varXBurn.slice(countBurn) = varX;
			paramXBurn.slice(countBurn) = paramX;
			paramPhyloBurn(countBurn) = paramPhylo;

			for(int j = 0; j < nRandom; j++){
				paramLatentBurn(countBurn,j) = paramLatent(j,0);
				latentBurn(countBurn,j) = latent(j,0);
			}

			// Counter
			countBurn++;
		}else if(i>=nburn && i%thin==0){
			// Save results for estimation
			meansParamXEst.col(countEst) = meansParamX;
			varXEst.slice(countEst) = varX;
			paramXEst.slice(countEst) = paramX;
			paramPhyloEst(countEst) = paramPhylo;

			for(int j = 0; j < nRandom; j++){
				paramLatentEst(countEst,j) = paramLatent(j,0);
				latentEst(countEst,j) = latent(j,0);
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
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXBurn)),
								   Rcpp::Named("meansParamX", wrap(trans(meansParamXBurn))),
							 	   Rcpp::Named("varX", wrap(varXBurn)),
							 	   Rcpp::Named("paramPhylo", wrap(paramPhyloBurn)),
							 	   Rcpp::Named("paramLatent", wrap(paramLatentBurn)),
							 	   Rcpp::Named("latent", wrap(latentBurn))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("meansParamX", wrap(trans(meansParamXEst))),
							 	   Rcpp::Named("varX", wrap(varXEst)),
							 	   Rcpp::Named("paramPhylo", wrap(paramPhyloEst)),
							 	   Rcpp::Named("paramLatent", wrap(paramLatentEst)),
							 	   Rcpp::Named("latent", wrap(latentEst))));

}
