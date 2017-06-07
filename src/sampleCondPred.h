#ifndef sampleCondPred_h
#define sampleCondPred_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::cube sampleCondPred(arma::mat& Y,
				   arma::mat& EstModel,
					 arma::vec residVar,
					 int nsite,
					 double nsp,
					 int nsample,
				 	 std::string family);

#endif
