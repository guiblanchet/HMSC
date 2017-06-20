#ifndef sampleCondPredAuto_h
#define sampleCondPredAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatentAuto.h"
#include "fixParamAuto.h"
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::cube sampleCondPredAuto(arma::mat& Y,
					 arma::field< arma::mat >& Auto,
					 arma::umat& RandomAuto,
					 arma::field< arma::mat >& latentAuto,
					 arma::field< arma::mat >& paramLatentAuto,
					 arma::field<arma::vec>& paramAuto,
					 arma::mat& residVar,
					 Rcpp::NumericMatrix& priorParamAutoDist,
					 int nsite,
					 double nsp,
					 int nAuto,
					 arma::vec& nAutoLev,
					 int npriorParamAuto,
					 int niter,
					 int nsample,
				 	 std::string family);
#endif
