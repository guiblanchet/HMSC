#ifndef rmvnorm_h
#define rmvnorm_h
#include <RcppArmadillo.h>
#include "rmvnorm.h"

arma::mat rmvnorm(int n, arma::vec Mean, arma::mat Var);

#endif
