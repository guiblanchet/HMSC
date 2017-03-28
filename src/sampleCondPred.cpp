#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleCondPred.h"

using namespace arma;
using namespace Rcpp;

// Calculates a prediction conditional on a subset of species.
//' @export
//[[Rcpp::export]]
RcppExport SEXP sampleCondPred(arma::mat& Y,
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
				YlatentTmp.slice(j) = sampleYlatentProbit(Y0Loc, Y1Loc, YlatentIni, EstModel.slice(j), residVar.row(j), nsp, nsite);
			}
			Ylatent(i,0) = YlatentTmp;
		}
	}

	if(family == "gaussian"){
		for (int i = 0; i < nsample; i++) {
			for(int j = 0; j < niter; j++){
				YlatentTmp.slice(i) = randn(nsite,nsp)*residVar+EstModel.slice(j);
			}
			Ylatent(i,0) = YlatentTmp;
		}
	}

	// return result object
	return Rcpp::List::create(Rcpp::Named("Ylatent", wrap(Ylatent)));
}
