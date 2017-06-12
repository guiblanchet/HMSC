#ifndef sampleCondPredLatentAuto_h
#define sampleCondPredLatentAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatent.h"
#include "updateLatentAuto.h"
#include "fixParamAuto.h"
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::field<arma::cube> sampleCondPredLatentAuto(arma::mat& Y,
					 arma::field< arma::mat >& Auto,
					 arma::umat& Random,
					 arma::umat& RandomAuto,
					 arma::field< arma::mat >& latent,
					 arma::field< arma::mat >& paramLatent,
					 arma::field< arma::mat >& latentAuto,
					 arma::field< arma::mat >& paramLatentAuto,
					 arma::field<arma::vec>& paramAuto,
					 arma::mat residVar,
					 Rcpp::NumericMatrix& priorParamAutoDist,
					 int nsite,
					 double nsp,
					 int nRandom,
					 arma::vec& nRandomLev,
					 int nAuto,
					 arma::vec& nAutoLev,
					 int npriorParamAuto,
					 int niter,
					 int nsample,
				 	 std::string family);
#endif
