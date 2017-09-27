#ifndef mcmcBinomialLogitX_h
#define mcmcBinomialLogitX_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentBinomialLogit.h"
#include "updateParamX.h"
#include "updatePrecX.h"
#include "updateMeansParamX.h"

RcppExport SEXP mcmcBinomialLogitX(arma::mat& Y,
							arma::mat& Ylatent,
							arma::mat& X,
							arma::mat& paramX,
							arma::mat& meansParamX,
							arma::mat& precX,
							arma::vec& residVar,
							arma::mat& priorMeansParamX,
							arma::mat& priorVarXScaleMat,
							double priorVarXDf,
							double priorResidVarScale,
							double priorResidVarShape,
							double ncount,
							double nsp,
							int nsite,
							int nparamX,
							int niter,
							int nburn,
							int thin,
							int verbose);

#endif
