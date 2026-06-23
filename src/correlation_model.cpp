// [[Rcpp::depends(BH, RcppParallel)]]
#include <Rcpp.h>
#include <cmath>

#include "rho_matern.h"

using namespace Rcpp;

// [[Rcpp::export]]
double matern(double d, double range, double smooth){
  return rho_matern(d, range, smooth);
}
