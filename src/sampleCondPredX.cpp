#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPredX.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//[[Rcpp::export]]
arma::field<arma::cube> sampleCondPredX(arma::mat& Y,
					 arma::mat& X,
					 arma::cube& paramX,
					 arma::vec residVar,
					 int nsite,
					 double nsp,
					 int niter,
					 int nsample,
				 	 std::string family) {

	// Define objects to store results
	cube YlatentSample(nsite, nsp, nsample);
	field<cube> Ylatent(niter,1);

	mat EstModel(nsite, nsp);

	for(int i = 0; i < niter ; i++){
		EstModel = X * trans(paramX.slice(i));

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
			// Store in final object
			Ylatent(j,0) = YlatentSample;
		}
	}
	// return result object
	return Ylatent;
}
