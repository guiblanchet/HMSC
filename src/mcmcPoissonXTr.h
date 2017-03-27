#ifndef mcmcPoissonXTr_h
#define mcmcPoissonXTr_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentPoisson.h"
#include "updateParamX.h"
#include "updatePrecXTr.h"
#include "updateParamTr.h"

RcppExport SEXP mcmcPoissonXTr(arma::mat& Y,
							  arma::mat& Ylatent,
							  arma::mat& X,
							  arma::mat& Tr,
							  arma::mat& paramX,
							  arma::mat& paramTr,
							  arma::mat& precX,
							  arma::vec& residVar,
							  arma::mat& priorParamTr,
							  arma::mat& priorVarTr,
							  arma::mat& priorVarXScaleMat,
							  double priorVarXDf,
							  double priorResidVarScale,
							  double priorResidVarShape,
							  double nsp,
							  int nsite,
							  int nparamX,
							  int nTr,
							  int niter,
							  int nburn,
							  int thin,
							  int verbose);

#endif
