#ifndef sampleYlatentPoisson_h
#define sampleYlatentPoisson_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm.h"

arma::mat sampleYlatentPoisson(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  arma::vec residVar,
							  double nsp,
							  int nsite);

#endif
