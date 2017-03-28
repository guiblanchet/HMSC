#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPred.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
arma::mat sampleCondPred(arma::mat& Y,
				   arma::cube& EstModel,
					 arma::vec residVar,
					 int niter,
					 double nsp,
					 int nsite,
				 	 std::string family) {

	// Define object to store results
	cube Ylatent(nsite, nsp, niter);
	if(family == "probit"){
		uvec Y0Loc = find(Y==0);
		uvec Y1Loc = find(Y==1);

		mat Ylatent = zeros(nsite,nsp);

		for (int i = 0; i < niter; i++) {
			Ylatent.slice(i) = sampleYlatentProbit(Y0Loc, Y1Loc, Ylatent, EstModel.slice(i), residVar, nsp, nsite)
		}
	}

	if(family == "gaussian"){
		for (int i = 0; i < niter; i++) {
			Ylatent.slice(i) = randn(nsite,nsp)*residVar+EstModel.slice(i)
		}
	}

	// return result object
	return Rcpp::Named("Ylatent", wrap(Ylatent));
}
