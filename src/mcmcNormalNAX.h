#ifndef mcmcNormalNAX_h
#define mcmcNormalNAX_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentNormal.h"
#include "updateParamX.h"
#include "updatePrecX.h"
#include "updateMeansParamX.h"
#include "updateResidVar.h"

RcppExport SEXP mcmcNormalNAX(arma::mat& Y,
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
							double nsp,
							int nsite,
							int nparamX,
							int niter,
							int nburn,
							int thin,
							int verbose);

#endif
