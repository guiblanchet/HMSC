#ifndef mcmcProbitAutoLatent_h
#define mcmcProbitAutoLatent_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"
#include "updateParamLatent.h"
#include "updateParamAuto.h"
#include "updateLatent.h"
#include "updateLatentAuto.h"
#include "updateShrinkLocal.h"
#include "updateParamShrinkGlobal.h"
#include "adaptVar.h"
#include "adaptVarAuto.h"
#include "fixParamAuto.h"

RcppExport SEXP mcmcProbitAutoLatent(arma::mat& Y,
									 arma::mat& Ylatent,
									 arma::field< arma::mat >& Auto,
									 arma::umat& RandomAuto,
									 arma::umat& Random,
									 arma::vec& residVar,
									 arma::field< arma::mat >& latent,
									 arma::field< arma::mat >& paramLatent,
									 arma::field< arma::mat >& shrinkLocal,
									 arma::field< arma::vec >& paramShrinkGlobal,
									 arma::field< arma::vec >& paramAuto,
									 arma::field< arma::mat >& latentAuto,
									 arma::field< arma::mat >& paramLatentAuto,
									 arma::field< arma::mat >& shrinkLocalAuto,
									 arma::field< arma::vec >& paramShrinkGlobalAuto,
									 double priorResidVarScale,
									 double priorResidVarShape,
									 arma::mat& priorParamAutoWeight,
									 Rcpp::NumericMatrix& priorParamAutoDist,
									 double priorShrinkLocal,
									 double priorShrinkOverallShape,
									 double priorShrinkOverallScale,
									 double priorShrinkSpeedShape,
									 double priorShrinkSpeedScale,
									 arma::vec& adapt,
									 arma::vec& redund,
									 int nAuto,
									 arma::vec& nAutoLev,
									 arma::vec& nLatentAuto,
									 int nRandom,
									 arma::vec& nRandomLev,
									 arma::vec& nLatent,
									 double nsp,
									 int nsite,
									 int npriorParamAuto,
									 int niter,
									 int nburn,
									 int thin,
									 int verbose);

#endif