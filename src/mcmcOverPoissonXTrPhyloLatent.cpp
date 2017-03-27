#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "mcmcOverPoissonXTrPhyloLatent.h"

using namespace arma ;
using namespace Rcpp ;

// Gibbs sampling for a model that includes a set of explanatory variable and latent variables.
//' @rdname mcmcProbitX
//' @export
// [[Rcpp::export]]
RcppExport SEXP mcmcOverPoissonXTrPhyloLatent(arma::mat& Y,
										 arma::mat& Ylatent,
										 arma::mat& X,
										 arma::mat& Tr,
										 arma::mat& Phylo,
										 arma::mat& iPhylo,
										 arma::umat& Random,
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
	cube paramTrBurn(nparamX,nTr,nburn/thin);
	cube paramXBurn(nsp,nparamX,nburn/thin);
	cube varXBurn(nparamX,nparamX,nburn/thin);

	vec paramPhyloBurn(nburn/thin);

	field<mat> latentBurn(nburn/thin,nRandom);
	field<mat> paramLatentBurn(nburn/thin,nRandom);

	mat varDistBurn(nsp,nburn/thin);

	// Define the result objects for estimation
	int nEst;
	nEst = niter-nburn;

	cube paramTrEst(nparamX,nTr,nEst/thin);
	cube paramXEst(nsp,nparamX,nEst/thin);
	cube varXEst(nparamX,nparamX,nEst/thin);

	vec paramPhyloEst(nEst/thin);

	field<mat> latentEst(nEst/thin,nRandom);
	field<mat> paramLatentEst(nEst/thin,nRandom);

	mat varDistEst(nsp,nEst/thin);

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

		// Recalculate meansParamX
		meansParamX = trans(paramTr*Tr);

		// Calculate the model estimate
		EstModel = X*trans(paramX);

		// Remove influence of X variables
		Yresid = Ylatent-EstModel;

		//------------------
		// Update paramLatent
		//------------------
		paramLatent = updateParamLatent(Yresid, Random, residVar, paramLatent, latent, shrink, nRandom, nLatent, nsp, nsite);

		//--------------
		// Update latent
		//--------------
		latent = updateLatent(Yresid, Random, residVar, paramLatent, latent, nRandom, nRandomLev, nLatent, nsp, nsite);

		//----------------
		// Sample Y latent
		//----------------
		// Include the effect of the latent variables
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

		//--------------
		// Update paramX
		//--------------
		// Remove the effect of the latent variables
		Yresid = Ylatent;
		for(int j = 0; j < nRandom; j++){
			mat latentMat = latent(j,0);
			Yresid = Yresid - latentMat.rows(Random.col(j))*trans(paramLatent(j,0));
		}

		paramX = updateParamXPhylo(Yresid, X, paramX, meansParamX, precX, residVar, wPhyloInvMat, nsp, nsite, nparamX);

		//--------------------------------
		// Update precX and calculate varX
		//--------------------------------
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

		//-------------
		// Save results
		//-------------
		if(i<nburn && i%thin==0){
			// Save burning results
			paramTrBurn.slice(countBurn) = paramTr;
			varXBurn.slice(countBurn) = precX.i();
			paramXBurn.slice(countBurn) = paramX;
			paramPhyloBurn(countBurn) = paramPhylo;
			for(int j = 0; j < nRandom; j++){
				paramLatentBurn(countBurn,j) = paramLatent(j,0);
				latentBurn(countBurn,j) = latent(j,0);
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
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXBurn)),
								   Rcpp::Named("paramTr", wrap(paramTrBurn)),
							 	   Rcpp::Named("varX", wrap(varXBurn)),
							 	   Rcpp::Named("paramPhylo", wrap(paramPhyloBurn)),
							 	   Rcpp::Named("paramLatent", wrap(paramLatentBurn)),
							 	   Rcpp::Named("latent", wrap(latentBurn)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistBurn)))),
				Rcpp::List::create(Rcpp::Named("paramX", wrap(paramXEst)),
								   Rcpp::Named("paramTr", wrap(paramTrEst)),
							 	   Rcpp::Named("varX", wrap(varXEst)),
							 	   Rcpp::Named("paramPhylo", wrap(paramPhyloEst)),
							 	   Rcpp::Named("paramLatent", wrap(paramLatentEst)),
							 	   Rcpp::Named("latent", wrap(latentEst)),
								   Rcpp::Named("varPoisson", wrap(trans(1/varDistEst)))));

}
