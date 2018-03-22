#ifndef mcmcBinomialLogitXPhylo_h
#define mcmcBinomialLogitXPhylo_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "sampleYlatentBinomialLogit.h"
#include "updateParamXPhylo.h"
#include "updateMeansParamXPhylo.h"
#include "updatePrecXPhylo.h"
#include "fixParamPhylo.h"
#include "updateParamPhylo.h"

RcppExport SEXP mcmcBinomialLogitXPhylo(arma::mat& Y,
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
								 double ncount,
								 int verbose);

#endif
