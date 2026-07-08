// [[Rcpp::depends(BH, RcppParallel)]]

#include <Rcpp.h>

#include <cmath>
#include <limits>

#include "bcl_storage.h"
#include "parallel_bcl.h"
#include "rho_matern.h"

using namespace Rcpp;


struct GaussianMarginalContribution {

  const std::vector<double> beta;

  double sill;
  double range;
  double t;
  double t2;
  double smooth;

  GaussianMarginalContribution(const Rcpp::NumericVector& beta_,
                               double sill_,
                               double range_,
                               double t_,
                               double t2_,
                               double smooth_)
    : beta(beta_.begin(), beta_.end()),
      sill(sill_),
      range(range_),
      t(t_),
      t2(t2_),
      smooth(smooth_) {

    if (beta.empty()) {
      Rcpp::stop("beta must have positive length.");
    }

    if (t_ <= 0.0 || !std::isfinite(t_)) {
      Rcpp::stop("t must be positive and finite.");
    }

    if (sill_ <= 0.0 || !std::isfinite(sill_)) {
      Rcpp::stop("sill must be positive and finite.");
    }

    if (range_ <= 0.0 || !std::isfinite(range_)) {
      Rcpp::stop("range must be positive and finite.");
    }

    if (smooth_ <= 0.0 || !std::isfinite(smooth_)) {
      Rcpp::stop("smooth must be positive and finite.");
    }

    if (t2_ <= 0.0 || !std::isfinite(t2_)) {
      Rcpp::stop("t2 must be positive and finite.");
    }
  }

  inline double operator()(double ei,
                         double ej,
                         double h) const {

    const double rho = rho_matern(h, range, smooth);
    const double r = rho * sill;
    const double s = t2 - r * r;

    if (s <= 0.0 || !std::isfinite(s)) {
      return 0.0;
    }

    return std::log(s) +
      (t * (ei * ei + ej * ej) - 2.0 * r * ei * ej) / s;
  }
};


struct GaussianConditionalContribution {

  const std::vector<double> beta;

  double sill;
  double range;
  double t;
  double t2;
  double smooth;
  double logt;

  GaussianConditionalContribution(const Rcpp::NumericVector& beta_,
                                  double sill_,
                                  double range_,
                                  double t_,
                                  double t2_,
                                  double smooth_)
    : beta(beta_.begin(), beta_.end()),
      sill(sill_),
      range(range_),
      t(t_),
      t2(t2_),
      smooth(smooth_),
      logt(0.0) {

    if (beta.empty()) {
      Rcpp::stop("beta must have positive length.");
    }

    if (t_ <= 0.0 || !std::isfinite(t_)) {
      Rcpp::stop("t must be positive and finite.");
    }

    if (sill_ <= 0.0 || !std::isfinite(sill_)) {
      Rcpp::stop("sill must be positive and finite.");
    }

    if (range_ <= 0.0 || !std::isfinite(range_)) {
      Rcpp::stop("range must be positive and finite.");
    }

    if (smooth_ <= 0.0 || !std::isfinite(smooth_)) {
      Rcpp::stop("smooth must be positive and finite.");
    }

    if (t2_ <= 0.0 || !std::isfinite(t2_)) {
      Rcpp::stop("t2 must be positive and finite.");
    }

    logt = std::log(t_);
  }

  inline double operator()(double ei,
                         double ej,
                         double h) const {

    const double rho = rho_matern(h, range, smooth);
    const double r = rho * sill;
    const double s = t2 - r * r;

    if (s <= 0.0 || !std::isfinite(s)) {
      return 0.0;
    }

    const double diff = ei - (r * ej / t);

    return std::log(s) - logt + (t / s) * diff * diff;
  }
};


// [[Rcpp::export]]
SEXP prepare_bcl_storage(Rcpp::List y_list_R,
                         Rcpp::List X_list_R,
                         Rcpp::List valid_pairs_list_R,
                         Rcpp::List dist_pairs_list_R) {

  BCLPtr ptr(
      new BCLStorage(
          y_list_R,
          X_list_R,
          valid_pairs_list_R,
          dist_pairs_list_R
      ),
      true
  );

  return ptr;
}


// [[Rcpp::export]]
double bcl_gaussian_marginal_ptr(SEXP ptr_,
                                 Rcpp::NumericVector beta,
                                 double sill,
                                 double range,
                                 double t,
                                 double t2,
                                 double smooth,
                                 int nthreads = 0,
                                 int min_pairs_parallel = 20000) {

  BCLPtr ptr(ptr_);

  GaussianMarginalContribution contrib(
      beta,
      sill,
      range,
      t,
      t2,
      smooth
  );

  const double out = 0.5 * parallel_block_pairwise_sum(
    *ptr,
    contrib,
    nthreads,
    min_pairs_parallel
  );

  if (!std::isfinite(out)) {
    return 1e20;
  }

  return out;
}


// [[Rcpp::export]]
double bcl_gaussian_conditional_ptr(SEXP ptr_,
                                    Rcpp::NumericVector beta,
                                    double sill,
                                    double range,
                                    double t,
                                    double t2,
                                    double smooth,
                                    int nthreads = 0,
                                    int min_pairs_parallel = 20000) {

  BCLPtr ptr(ptr_);

  GaussianConditionalContribution contrib(
      beta,
      sill,
      range,
      t,
      t2,
      smooth
  );

  const double out = 0.5 * parallel_block_pairwise_sum(
    *ptr,
    contrib,
    nthreads,
    min_pairs_parallel
  );

  if (!std::isfinite(out)) {
    return 1e20;
  }

  return out;
}
