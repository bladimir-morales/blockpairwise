// [[Rcpp::depends(BH, RcppParallel)]]
#include <Rcpp.h>
#include <cmath>

#include "rho_matern.h"

using namespace Rcpp;

// [[Rcpp::export]]
double matern(double d, double phi, double nu){
  return rho_matern(d, phi, nu);
}
