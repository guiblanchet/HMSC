#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPredLatent.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//[[Rcpp::export]]
arma::cube sampleCondPredLatent(arma::mat& Y,
					 arma::umat& Random,
					 arma::field< arma::mat >& latent,
					 arma::field< arma::mat >& paramLatent,
					 arma::mat residVar,
					 int nsite,
					 double nsp,
					 int nRandom,
					 arma::vec& nRandomLev,
					 int niter,
					 int nsample,
				 	 std::string family) {

	// Define basic objects to store results
	mat YlatentSample(nsite, nsp);
	cube Ylatent(nsite, nsp, niter);
	mat EstModel(nsite, nsp);

	mat residVarT = trans(residVar);

	// Define field object to store one iterations of latent variables and their associated parameters
	field<mat> latent1iter(nRandom,1);
	field<mat> paramLatent1iter(nRandom,1);
	vec nLatent(nRandom);

	// Reorganize latent
	field<mat> latentOrg(niter, nRandom);
	field<mat> paramLatentOrg(niter, nRandom);
	int counter = 0;
	for(int j = 0; j < nRandom ; j++){
		for(int i = 0; i < niter ; i++){
			latentOrg(i,j) = latent(counter,0);
			paramLatentOrg(i,j) = paramLatent(counter,0);
			counter++;
		}
	}

	for(int i = 0; i < niter ; i++){
		// Calculate the model estimation
		EstModel.zeros();

		for(int j = 0; j < nRandom ; j++){
			mat latentMat = latentOrg(i,j);
			nLatent(j) = latentMat.n_cols;
			mat paramLatentMat = paramLatentOrg(i,j);
			EstModel = EstModel + latentMat.rows(Random.col(j))*trans(paramLatentMat);
		}

		// Sample Ylatent
		for (int k = 0; k < nsample; k++) {
			if(family == "probit"){
				uvec Y0Loc = find(Y==0);
				uvec Y1Loc = find(Y==1);
				uvec YNALoc = find_nonfinite(Y);

				mat YlatentIni = zeros(nsite,nsp);

				YlatentSample = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, YlatentIni, EstModel, residVarT.col(i), nsp, nsite);
			}

			if(family == "gaussian"){
				mat repResidVar(nsite,nsp);
				repResidVar = repmat(residVarT.col(i),nsite,1);
				YlatentSample = randn(nsite,nsp)%repResidVar+EstModel;
			}

			if(family == "poisson" | family == "overPoisson"){
				mat YlatentIni = zeros(nsite,nsp);
				YlatentSample = sampleYlatentPoisson(Y, YlatentIni, EstModel, residVarT.col(i), nsp, nsite);
			}

			// Update latent
			for(int j = 0; j < nRandom ; j++){
				latent1iter(j,0) = latentOrg(i,j);
				paramLatent1iter(j,0) = paramLatentOrg(i,j);
			}

			latent1iter = updateLatent(YlatentSample,Random,residVarT.col(i),paramLatent1iter,latent1iter,nRandom,nRandomLev,nLatent,nsp,nsite);

			for(int j = 0; j < nRandom ; j++){
				latentOrg(i,j) = latent1iter(j,0);
			}
		}
		// Save results
		Ylatent.slice(i) = YlatentSample;
	}

	// return result object
	return Ylatent;
}
