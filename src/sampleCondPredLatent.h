#ifndef sampleCondPredLatent_h
#define sampleCondPredLatent_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatent.h"
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::field<arma::cube> sampleCondPredLatent(arma::mat& Y,
					 arma::umat& Random,
					 arma::field< arma::mat >& latent,
					 arma::field< arma::mat >& paramLatent,
					 arma::vec residVar,
					 int nsite,
					 double nsp,
					 int nRandom,
					 arma::vec& nRandomLev,
					 int niter,
					 int nsample,
				 	 std::string family);
#endif
