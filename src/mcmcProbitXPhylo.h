#ifndef mcmcProbitXPhylo_h
#define mcmcProbitXPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentProbit.h"
#include "updateParamXPhylo.h"
#include "updateMeansParamXPhylo.h"
#include "updatePrecXPhylo.h"
#include "fixParamPhylo.h"
#include "updateParamPhylo.h"

RcppExport SEXP mcmcProbitXPhylo(arma::mat& Y,
								 arma::mat& Ylatent,
								 arma::mat& X,
								 arma::mat& Phylo,
								 arma::mat& iPhylo,
								 arma::mat& paramX,
								 arma::mat& meansParamX,
								 double paramPhylo,
								 arma::mat& precX,
								 arma::vec& residVar,
								 arma::mat& priorMeansParamX,
								 arma::mat& priorVarMeansParamX,
								 arma::mat& priorVarXScaleMat,
								 double priorVarXDf,
								 double priorResidVarScale,
								 double priorResidVarShape,
								 arma::mat& priorParamPhyloWeight,
								 Rcpp::NumericVector& priorParamPhyloGrid,
								 double nsp,
								 int nsite,
								 int nparamX,
								 int nparamPhylo,
								 int niter,
								 int nburn,
								 int thin,
								 int verbose);

#endif
