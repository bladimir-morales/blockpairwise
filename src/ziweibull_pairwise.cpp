// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH, RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>

#include <cmath>
#include <algorithm>
#include <vector>
#include <thread>

#include "bcl_zi_storage.h"
#include "parallel_bcl_zi.h"
#include "rho_matern.h"
#include <R_ext/RS.h>

extern "C" {
  void F77_CALL(zibvn)(
      const double* a1,
      const double* a2,
      const double* rho,
      double* p
  );
}

using namespace Rcpp;

// ============================================================
// Bivariate zero-zero probability
// ============================================================

inline double log_p00_cpp(double a1, double a2, double rho) {

  if (!std::isfinite(a1) || !std::isfinite(a2) ||
      !std::isfinite(rho) || std::abs(rho) >= 1.0) {
      return NAN;
  }

  double p00 = NAN;

  F77_CALL(zibvn)(&a1, &a2, &rho, &p00);

  if (!std::isfinite(p00) || p00 < 0.0 || p00 > 1.0) {
    return NAN;
  }

  if (p00 == 0.0) return -INFINITY;

  return std::log(p00);
}

// ============================================================
// Log standard normal CDF
// ============================================================

inline double log_pnorm_zi_cpp(double x) {

  if (std::isnan(x)) return NAN;
  if (x == -INFINITY) return -INFINITY;
  if (x == INFINITY) return 0.0;

  const double inv_sqrt2 = 0.70710678118654752440;

  if (x >= 0.0) {
    return std::log1p(-0.5 * std::erfc(x * inv_sqrt2));
  }

  if (x > -10.0) {
    return std::log(0.5 * std::erfc(-x * inv_sqrt2));
  }

  // Mills ratio: continued-fraction evaluation in the lower tail
  const double t = -x;
  double cf = 0.0;

  for (int k = 40; k >= 1; --k) {
    cf = static_cast<double>(k) / (t + cf);
  }

  return -0.5 * t * t - 0.91893853320467274178 -
    std::log(t + cf);
}


// ============================================================
// Gaussian copula log-density
// ============================================================

inline double log_c_gauss_zi_cpp(double b1, double b2, double rho) {

  if (!std::isfinite(b1) || !std::isfinite(b2) ||
      !std::isfinite(rho) || std::abs(rho) >= 1.0) {
      return NAN;
  }

  const double den = (1.0 - rho) * (1.0 + rho);

  return -0.5 * std::log(den) +
    (rho * b1 * b2 -
    0.5 * rho * rho * (b1 * b1 + b2 * b2)) / den;
}

// ============================================================
// Zero-inflated contribution with precomputed marginals
// ============================================================

struct BCLZIMarginalContribution {

  double range;
  double smooth;

  BCLZIMarginalContribution(double range_, double smooth_)
    : range(range_), smooth(smooth_) {

    if (!std::isfinite(range) || range <= 0.0 ||
        !std::isfinite(smooth) || smooth <= 0.0) {
        Rcpp::stop("Invalid covariance parameters.");
    }
  }

  inline double operator()(double y1, double y2,
                         double a1, double a2,
                         double b1, double b2,
                         double L1, double L2,
                         double d) const {

    const double rho = rho_matern(d, range, smooth);

    if (!std::isfinite(rho) || std::abs(rho) >= 1.0) {
      return NAN;
    }

    const double den = (1.0 - rho) * (1.0 + rho);
    const double sd = std::sqrt(den);

    // Case 00
    if (y1 == 0.0 && y2 == 0.0) {
      return log_p00_cpp(a1, a2, rho);
    }

    // Case 0+
    if (y1 == 0.0 && y2 > 0.0) {
      return L2 + log_pnorm_zi_cpp((a1 - rho * b2) / sd);
    }

    // Case +0
    if (y1 > 0.0 && y2 == 0.0) {
      return L1 + log_pnorm_zi_cpp((a2 - rho * b1) / sd);
    }

    // Case ++
    if (y1 > 0.0 && y2 > 0.0) {
      return L1 + L2 + log_c_gauss_zi_cpp(b1, b2, rho);
    }

    return NAN;
  }
};




// ============================================================
// Exported storage constructor
// ============================================================

// [[Rcpp::export]]
SEXP prepare_zi_storage_cpp(List y_list,
                            List pairs_list,
                            List dist_list) {

  BCLZIStoragePtr ptr(
      new BCLZIStorage(
          y_list,
          pairs_list,
          dist_list
      ),
      true
  );

  return ptr;
}


// ============================================================
// Objective with precomputed ZI marginals
// ============================================================

// [[Rcpp::export]]
double bp_zi_marginals_ptr(SEXP ptr_,
                           Rcpp::List marg_blocks,
                           double range,
                           double smooth,
                           int nthreads = 0) {

  BCLZIStoragePtr ptr(ptr_);

  BCLZIMarginals marg(marg_blocks, *ptr);

  BCLZIMarginalContribution contrib(range, smooth);

  const double cl = parallel_zi_marginals_block_sum(
    *ptr, marg, contrib, nthreads
  );

  return -cl;
}

