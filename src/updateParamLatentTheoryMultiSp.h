#ifndef updateParamLatentTheoryMultiSp_h
#define updateParamLatentTheoryMultiSp_h
#include <RcppArmadillo.h>
#include <Rcpp.h>

arma::field<arma::mat> updateParamLatentTheoryMultiSp(arma::mat& Ylatent,
										 arma::umat& Random,
										 arma::vec& residVar,
										 arma::field<arma::mat>& paramLatent,
										 arma::field<arma::mat>& latent,
										 arma::field<arma::mat>& shrink,
										 int nRandom,
										 arma::vec& nLatent,
										 double nsp,
										 int nsite,
									 	 arma::mat& diagMat);

#endif
