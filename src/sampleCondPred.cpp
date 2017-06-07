#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPred.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//[[Rcpp::export]]
arma::cube sampleCondPred(arma::mat& Y,
				   arma::mat& EstModel,
					 arma::vec residVar,
					 int nsite,
					 double nsp,
					 int nsample,
				 	 std::string family) {

	// Define objects to store results
//	field<cube> Ylatent(nsample,1);
	cube Ylatent(nsite, nsp, nsample);

	if(family == "probit"){
		uvec Y0Loc = find(Y==0);
		uvec Y1Loc = find(Y==1);
		uvec YNALoc = find_nonfinite(Y);

		mat YlatentIni = zeros(nsite,nsp);

		for (int i = 0; i < nsample; i++) {
					Ylatent.slice(i) = sampleYlatentProbit(Y0Loc, Y1Loc, YNALoc, YlatentIni, EstModel, residVar, nsp, nsite);
		}
	}

	if(family == "gaussian"){
		mat repResidVar(nsite,nsp);
		for (int i = 0; i < nsample; i++) {
				repResidVar = repmat(residVar,nsite,1);
				Ylatent.slice(i) = randn(nsite,nsp)%repResidVar+EstModel;
		}
	}

	if(family == "poisson" | family == "overPoisson"){
		mat YlatentIni = zeros(nsite,nsp);

		for (int i = 0; i < nsample; i++) {
				Ylatent.slice(i) = sampleYlatentPoisson(Y, YlatentIni, EstModel, residVar, nsp, nsite);
		}
	}

	// return result object
	return Ylatent;
}
