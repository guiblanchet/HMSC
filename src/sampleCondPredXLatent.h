#ifndef sampleCondPredXLatent_h
#define sampleCondPredXLatent_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatent.h"
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::cube sampleCondPredXLatent(arma::mat& Y,
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
				 	 std::string family);
#endif
