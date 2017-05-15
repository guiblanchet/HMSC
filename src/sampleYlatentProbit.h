#ifndef sampleYlatentProbit_h
#define sampleYlatentProbit_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm.h"

arma::mat sampleYlatentProbit(arma::uvec& Y0Loc,
								arma::uvec& Y1Loc,
								arma::uvec& YNALoc,
							  arma::mat& Ylatent,
							  arma::mat& EstModel,
							  arma::vec& residVar,
							  double nsp,
							  int nsite);

#endif
