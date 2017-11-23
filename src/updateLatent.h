#ifndef updateLatent_h
#define updateLatent_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::field<arma::mat> updateLatent(arma::mat& Ylatent,
									arma::umat& Random,
									arma::vec& residVar,
									arma::field<arma::mat>& paramLatent,
									arma::field<arma::mat>& latent,
									int nRandom,
									arma::vec& nRandomLev,
									arma::vec& nLatent,
									double nsp,
									int nsite);

#endif
