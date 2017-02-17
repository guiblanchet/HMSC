#ifndef mcmcProbitXLatent_h
#define mcmcProbitXLatent_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"
#include "updateParamX.h"
#include "updatePrecX.h"
#include "updateMeansParamX.h"
#include "updateParamLatent.h"
#include "updateLatent.h"
#include "updateShrinkLocal.h"
#include "updateParamShrinkGlobal.h"
#include "adaptVar.h"

RcppExport SEXP mcmcProbitXLatent(arma::mat& Y,
								  arma::mat& Ylatent,
								  arma::mat& X,
								  arma::umat& Random, 
								  arma::mat& paramX,
								  arma::mat& meansParamX,
								  arma::mat& precX,
								  arma::vec& residVar,
								  arma::mat& priorMeansParamX,
								  arma::mat& priorVarXScaleMat,
								  double priorVarXDf,
								  double priorResidVarScale,
								  double priorResidVarShape,
								  double priorShrinkLocal,
								  double priorShrinkOverallShape,
								  double priorShrinkOverallScale,
								  double priorShrinkSpeedShape,
								  double priorShrinkSpeedScale,
								  arma::vec& adapt,
								  arma::vec& redund,
								  int nRandom,
								  arma::vec& nRandomLev,
								  arma::vec& nLatent,
								  double nsp,
								  int nsite,
								  int nparamX,
								  int niter,
								  int nburn,
								  int thin,
								  int verbose);
							
#endif
