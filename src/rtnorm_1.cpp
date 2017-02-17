#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm_1.h"

using namespace arma ;
using namespace Rcpp ;

// Simplified sample from a truncated normal distribution where the standard deviation is always 1.
double rtnorm_1(double mean, double a, double b){
    // Initialize basic objects
    double lo = 0.0;
    double hi = 0.0;
    double out;
    double zlo = a-mean;
    double zhi = b-mean;
    
    if(a<b) {
			if ((std::fabs(zlo) < 8.0) && (std::fabs(zhi) < 8.0)) {
				lo = Rcpp::stats::pnorm_1(a, mean, true, false);
				hi = Rcpp::stats::pnorm_1(b, mean, true, false);
			} else if (std::fabs(zlo) >8.0) {
				lo = 0.0;
				hi = Rcpp::stats::pnorm_1(b, mean, true, false);
			} else if (std::fabs(zhi) >8.0) {
				lo = Rcpp::stats::pnorm_1(a, mean, true, false);
				hi = 1.0;
			}
			
			double u = Rf_runif(lo, hi);
			u = std::max(u, 6.7e-16);
			u = std::min(u, 1.0-6.7e-16);
	
			out = Rcpp::stats::qnorm_1(u, mean, true, false);
			//if inverse cdf fails numerically, fall back to slice sampler
			if (!arma::is_finite(out) || out<a || out>b) {
				out = rtnorm_slice(20, mean, a, b);
			}
    } else {
    	out=a;
    }
    return(out);
}
