#ifndef sampleYlatentBinomialLogit_h
#define sampleYlatentBinomialLogit_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm.h"

arma::mat sampleYlatentBinomialLogit(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  double ncount,
							  double nsp,
							  int nsite);

#endif
