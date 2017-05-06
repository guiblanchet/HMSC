#ifndef updateLatentAutoTheory1Sp_h
#define updateLatentAutoTheory1Sp_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rmvnorm.h"

arma::field<arma::mat> updateLatentAutoTheory1Sp(arma::mat& Yresid,
										arma::umat& RandomAuto,
										arma::vec& residVar,
										arma::field<arma::vec>& paramAuto,
										arma::field<arma::cube>& wAutoInv,
										arma::field<arma::mat>& paramLatentAuto,
										arma::field<arma::mat>& latentAuto,
										arma::mat& priorParamAutoDistArma,
										int nAuto,
										arma::vec& nAutoLev,
										arma::vec& nLatentAuto,
										double nsp,
										int nsite,
										arma::mat& diagMat);

#endif
