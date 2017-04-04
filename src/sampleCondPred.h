#ifndef sampleCondPred_h
#define sampleCondPred_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::field<arma::cube> sampleCondPred(arma::mat& Y,
				   arma::cube& EstModel,
					 arma::mat residVar,
					 int nsite,
					 double nsp,
					 int niter,
					 int nsample,
				 	 std::string family);

#endif
