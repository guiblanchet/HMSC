#ifndef updateLatentAuto_h
#define updateLatentAuto_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::field<arma::mat> updateLatentAuto(arma::mat& Yresid,
										arma::umat& RandomAuto,
										arma::vec residVar,
										arma::field<arma::vec>& paramAuto,
										arma::field<arma::cube>& wAutoInv,
										arma::field<arma::mat>& paramLatentAuto,
										arma::field<arma::mat>& latentAuto,
										arma::mat& priorParamAutoDistArma,
										int nAuto,
										arma::vec& nAutoLev,
										arma::vec nLatentAuto,
										double nsp,
										int nsite);

#endif
