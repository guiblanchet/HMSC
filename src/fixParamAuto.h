#ifndef fixParamAuto_h
#define fixParamAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "fixParamAuto.h"
#include "euclDist.h"

arma::field<arma::cube> fixParamAuto(arma::mat& Auto,
									 Rcpp::NumericMatrix& priorParamAutoDist,
									 double nsp, 
									 arma::vec& nAutoLev,
									 int npriorParamAuto,
									 int i);

#endif