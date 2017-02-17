#ifndef updateParamAuto_h
#define updateParamAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleNum.h"

// This function updates paramPhylo when traits are considered
arma::field<arma::vec> updateParamAuto(arma::field<arma::mat>& latentAuto,
									   arma::field<arma::vec>& paramAuto,
									   int npriorParamAuto,
									   arma::mat wAutoDet,
									   arma::field< arma::cube >& wAutoInv,
									   arma::mat& priorParamAutoWeight,
									   Rcpp::NumericMatrix& priorParamAutoDist,
									   arma::vec& nLatentAuto,
									   int nAuto);

#endif
