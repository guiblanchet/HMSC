#ifndef rtnorm_h
#define rtnorm_h
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm.h"

double rtnorm(double mean,
				double sd,
				double low,
				double high);
#endif
