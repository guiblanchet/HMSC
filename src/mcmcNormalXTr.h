#ifndef mcmcNormalXTr_h
#define mcmcNormalXTr_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "updateParamX1.h"
#include "updatePrecXTr.h"
#include "updateParamTr.h"
#include "updateResidVar.h"

RcppExport SEXP mcmcNormalXTr(arma::mat& Ylatent,
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
