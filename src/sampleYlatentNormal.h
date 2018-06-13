#ifndef sampleYlatentNormal_h
#define sampleYlatentNormal_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::mat sampleYlatentNormal(arma::mat& Y,
							  arma::mat& EstModel,
							  arma::vec& residVar,
							  double nsp,
							  int nsite);

#endif
