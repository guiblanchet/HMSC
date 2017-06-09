#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPredAuto.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//[[Rcpp::export]]
arma::field<arma::cube> sampleCondPredAuto(arma::mat& Y,
					 arma::field< arma::mat >& Auto,
					 arma::umat& RandomAuto,
					 arma::field< arma::mat >& latentAuto,
					 arma::field< arma::mat >& paramLatentAuto,
					 arma::field<arma::vec>& paramAuto,
					 arma::vec residVar,
					 Rcpp::NumericMatrix& priorParamAutoDist,
					 int nsite,
					 double nsp,
					 int nAuto,
					 arma::vec& nAutoLev,
					 int npriorParamAuto,
					 int niter,
					 int nsample,
				 	 std::string family) {

	// Define basic objects to store results
	cube YlatentSample(nsite, nsp, nsample);
	field<cube> Ylatent(nsite, nsp, nsample);
	mat EstModel(nsite, nsp);
	EstModel.zeros();

	// Define field object to store one iterations of latent variables and their associated parameters
	field<mat> latentAuto1iter(nAuto,1);
	field<mat> paramLatentAuto1iter(nAuto,1);
	vec nLatentAuto(nAuto);

	//Define wAutoInv object
	field<cube> wAutoInv(nAuto,1);
	field<cube> paramAutoNoGibb(2,1);
	for (int i = 0; i < nAuto; i++) {
		paramAutoNoGibb = fixParamAuto(Auto(i,0),priorParamAutoDist,nsp,nAutoLev,npriorParamAuto,i);
		wAutoInv(i,0) = paramAutoNoGibb(1,0);
	}

	// Convert NumericMatrix to arma::mat
	mat priorParamAutoDistArma = as<arma::mat>(priorParamAutoDist);

	////////////////////////////////
	// Sample conditional prediction
	////////////////////////////////
	for(int i = 0; i < niter ; i++){
		// Calculate the model estimation
		EstModel.zeros();

		for(int j = 0; j < nAuto ; j++){
			mat latentAutoMat = latentAuto(i,j);
			nLatentAuto(j) = latentAutoMat.n_cols;
			mat paramLatentAutoMat = paramLatentAuto(i,j);
			EstModel = EstModel + latentAutoMat.rows(RandomAuto.col(j))*trans(paramLatentAutoMat);
		}

		/////////////////
		// Sample Ylatent
		/////////////////
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

			////////////////
			// Update latent
			////////////////
			for(int j = 0; j < nAuto ; j++){
				latentAuto1iter(j,0) = latentAuto(i,j);
				paramLatentAuto1iter(j,0) = paramLatentAuto(i,j);
			}

			latentAuto1iter = updateLatentAuto(YlatentSample.slice(j), RandomAuto, residVar, paramAuto, wAutoInv, paramLatentAuto1iter, latentAuto1iter, priorParamAutoDistArma, nAuto, nAutoLev, nLatentAuto, nsp, nsite);

			for(int j = 0; j < nAuto ; j++){
				latentAuto(i,j) = latentAuto1iter(j,0);
			}
		}
	}

	// return result object
	return Ylatent;
}
