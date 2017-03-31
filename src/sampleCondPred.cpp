#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPred.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//' @export
//[[Rcpp::export]]
arma::field<arma::cube> sampleCondPred(arma::mat& Y,
				   arma::cube& EstModel,
					 arma::mat residVar,
					 int nsite,
					 double nsp,
					 int niter,
					 int nsample,
				 	 std::string family) {

	// Define objects to store results
	field<cube> Ylatent(nsample,1);
	cube YlatentTmp(nsite, nsp, niter);

	if(family == "probit"){
		uvec Y0Loc = find(Y==0);
		uvec Y1Loc = find(Y==1);

		mat YlatentIni = zeros(nsite,nsp);

		for (int i = 0; i < nsample; i++) {
			for(int j = 0; j < niter; j++){
				YlatentTmp.slice(j) = sampleYlatentProbit(Y0Loc, Y1Loc, YlatentIni, EstModel.slice(j), trans(residVar.row(j)), nsp, nsite);
			}
			Ylatent(i,0) = YlatentTmp;
		}
	}

	if(family == "gaussian"){
		mat repResidVar(nsite,nsp);
		for (int i = 0; i < nsample; i++) {
			for(int j = 0; j < niter; j++){
				repResidVar = repmat(residVar.row(j),nsite,1);
				YlatentTmp.slice(j) = randn(nsite,nsp)%repResidVar+EstModel.slice(j);
			}
			Ylatent(i,0) = YlatentTmp;
		}
	}

	// return result object
	return Ylatent;
}
