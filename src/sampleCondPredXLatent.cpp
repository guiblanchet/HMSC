#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPredXLatent.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//' @export
//[[Rcpp::export]]
arma::field<arma::cube> sampleCondPredXLatent(arma::mat& Y,
					 arma::mat& X,
					 arma::umat& Random,
					 arma::cube& paramX,
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
	cube YlatentSample(nsite, nsp, niter);
	field<cube> Ylatent(niter,1);
	mat EstModel(nsite, nsp);

	mat Yresid(nsite,nsp);

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
		EstModel = X * trans(paramX.slice(i));

		for(int j = 0; j < nRandom ; j++){
			mat latentMat = latentOrg(i,j);
			nLatent(j) = latentMat.n_cols;
			mat paramLatentMat = paramLatentOrg(i,j);
			EstModel = EstModel + latentMat.rows(Random.col(j))*trans(paramLatentMat);
		}

		// Sample Ylatent
		for (int j = 0; j < nsample; j++) {
			if(family == "probit"){
				uvec Y0Loc = find(Y==0);
				uvec Y1Loc = find(Y==1);
				uvec YNALoc = find_nonfinite(Y);

				mat YlatentIni = zeros(nsite,nsp);

				YlatentSample.slice(j) = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, YlatentIni, EstModel, residVarT.col(i), nsp, nsite);
			}

			if(family == "gaussian"){
				mat repResidVar(nsite,nsp);
				repResidVar = repmat(residVar.row(i),nsite,1);
				YlatentSample.slice(j) = randn(nsite,nsp)%repResidVar+EstModel;
			}

			if(family == "poisson" | family == "overPoisson"){
				mat YlatentIni = zeros(nsite,nsp);
				YlatentSample.slice(j) = sampleYlatentPoisson(Y, YlatentIni, EstModel, residVarT.col(i), nsp, nsite);
			}

			// Update latent
			for(int j = 0; j < nRandom ; j++){
				latent1iter(j,0) = latentOrg(i,j);
				paramLatent1iter(j,0) = paramLatentOrg(i,j);
			}

			// Remove influence of X variables
			Yresid = YlatentSample.slice(j) - X * trans(paramX.slice(i));

			latent1iter = updateLatent(Yresid,Random,residVarT.col(i),paramLatent1iter,latent1iter,nRandom,nRandomLev,nLatent,nsp,nsite);

			for(int j = 0; j < nRandom ; j++){
				latentOrg(i,j) = latent1iter(j,0);
			}
		}
		// Save results
		Ylatent(i,0) = YlatentSample;
	}

	// return result object
	return Ylatent;
}
