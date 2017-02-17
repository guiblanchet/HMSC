#include <RcppArmadillo.h>
#include <Rcpp.h>
#include "rtnorm_slice.h"

using namespace arma ;
using namespace Rcpp ;

double rtnorm_slice(int iter, double mean, double a, double b) {
	
	double up  = 0.0;
	double upp = 0.0;
	double lo = 0.0;
	double loo = 0.0;
	double y  = 0.0;
	double x = 0.0;
	
	// What to do with infinite value
	if(b==R_PosInf){
		if(a==-R_PosInf){
			x = 0;
		}else{
			x = a+Rf_runif(0, 1.0);
		}
	}else{
		if(a==-R_PosInf){
			x = b-Rf_runif(0, 1.0);
		}else{
			x = (b-a)/2;
		}
	}
	/////////////////
	// slice sampling
	/////////////////
	// When a is any value (except infinite)
	if(a!=-R_PosInf){
		for (int i=0; i<iter; i++) {
			y = exp(-0.5*(x-mean)*(x-mean))*Rf_runif(0.0, 1.0);
			upp = mean + sqrt(-2.0*log(y));
			up = (b==R_PosInf || upp<b) ? upp : b;
			lo = std::max(mean - sqrt(-2.0*log(y)), a);
			x = (up-lo)*Rf_runif(0.0, 1.0) + lo;
		}
	// When a is infinite
	}else{
		for (int i=0; i<iter; i++) {
			y = exp(-0.5*(x-mean)*(x-mean))*Rf_runif(0.0, 1.0);
			up = std::min(mean + sqrt(-2.0*log(y)), b);
			loo = mean - sqrt(-2.0*log(y));
			lo = (a==-R_PosInf || loo>a) ? loo : a;
			x = (lo-up)*Rf_runif(0.0, 1.0) + up;
		}
	}
	
	return(x);
	
}