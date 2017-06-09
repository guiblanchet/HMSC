#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPredLatent.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//[[Rcpp::export]]
arma::field<arma::cube> sampleCondPredLatent(arma::mat& Y,
					 arma::umat& Random,
					 arma::field< arma::mat >& latent,
					 arma::field< arma::mat >& paramLatent,
					 arma::vec residVar,
					 int nsite,
					 double nsp,
					 int nRandom,
					 arma::vec& nRandomLev,
					 int niter,
					 int nsample,
				 	 std::string family) {

	// Define basic objects to store results
	cube YlatentSample(nsite, nsp, nsample);
	field<cube> Ylatent(nsite, nsp, nsample);
	mat EstModel(nsite, nsp);

	// Define field object to store one iterations of latent variables and their associated parameters
	field<mat> latent1iter(nRandom,1);
	field<mat> paramLatent1iter(nRandom,1);
	vec nLatent(nRandom);

	for(int i = 0; i < niter ; i++){
		// Calculate the model estimation
		EstModel.zeros();

		for(int j = 0; j < nRandom ; j++){
			mat latentMat = latent(i,j);
			nLatent(j) = latentMat.n_cols;
			mat paramLatentMat = paramLatent(i,j);
			EstModel = EstModel + latentMat.rows(Random.col(j))*trans(paramLatentMat);
		}

		// Sample Ylatent
		for (int j = 0; j < nsample; j++) {
			if(family == "probit"){
				uvec Y0Loc = find(Y==0);
				uvec Y1Loc = find(Y==1);
				uvec YNALoc = find_nonfinite(Y);

				mat YlatentIni = zeros(nsite,nsp);

				YlatentSample.slice(j) = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, YlatentIni, EstModel, residVar, nsp, nsite);
			}

			if(family == "gaussian"){
				mat repResidVar(nsite,nsp);
				repResidVar = repmat(residVar,nsite,1);
				YlatentSample.slice(j) = randn(nsite,nsp)%repResidVar+EstModel;
			}

			if(family == "poisson" | family == "overPoisson"){
				mat YlatentIni = zeros(nsite,nsp);
				YlatentSample.slice(j) = sampleYlatentPoisson(Y, YlatentIni, EstModel, residVar, nsp, nsite);
			}

			// Update latent
			for(int j = 0; j < nRandom ; j++){
				latent1iter(j,0) = latent(i,j);
				paramLatent1iter(j,0) = paramLatent(i,j);
			}

			latent1iter = updateLatent(YlatentSample.slice(j),Random,residVar,paramLatent1iter,latent1iter,nRandom,nRandomLev,nLatent,nsp,nsite);

			for(int j = 0; j < nRandom ; j++){
				latent(i,j) = latent1iter(j,0);
			}
		}
	}

	// return result object
	return Ylatent;
}
