#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm.h"

using namespace arma ;
using namespace Rcpp ;

// Sample from a truncated normal distribution.

double rtnorm(double mean,
				double sd,
				double low,
				double high){
	
	// Sample from a uniform distribution
	double unifSample = Rf_runif(0, 1);
	
	// Define boundaries
	double lower = 2*Rcpp::stats::pnorm_0((low-mean)/sd,true,false)-1; // This is the same as : erf((low-mean)/(sqrt(2)*sd))
	double higher = 2*Rcpp::stats::pnorm_0((-mean+high)/sd,true,false)-1; // This is the same as : erf((-mean+high)/(sqrt(2)*sd))
	
	// Sample from a truncated normal distribution
	double sample = mean+sd*Rcpp::stats::qnorm_0((1+(1-unifSample)*lower+unifSample*higher)/2,true,false); // This is the same as : erfinv((1-unifSample)*lower+unifSample*higher)
	
	return(sample);
}
