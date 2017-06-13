#ifndef sampleCondPredXLatentAuto_h
#define sampleCondPredXLatentAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateLatent.h"
#include "updateLatentAuto.h"
#include "fixParamAuto.h"
#include "sampleYlatentProbit.h"
#include "sampleYlatentPoisson.h"

arma::cube sampleCondPredXLatentAuto(arma::mat& Y,
					 arma::mat& X,
					 arma::field< arma::mat >& Auto,
					 arma::umat& Random,
					 arma::umat& RandomAuto,
					 arma::cube& paramX,
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
					 int nAuto,
					 arma::vec& nAutoLev,
					 int npriorParamAuto,
					 int niter,
					 int nsample,
				 	 std::string family);
#endif
