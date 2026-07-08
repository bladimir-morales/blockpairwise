// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH, RcppParallel)]]

#include <Rcpp.h>
#include <RcppParallel.h>

#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <boost/math/distributions/normal.hpp>

#include <cmath>
#include <algorithm>
#include <vector>
#include <thread>

#include "bcl_zi_storage.h"
#include "parallel_bcl_zi.h"
#include "rho_matern.h"

using namespace Rcpp;


// ============================================================
// Pure C++ utilities
// ============================================================

inline double clamp01_cpp(double u,
                          const double eps = 1e-10) {

  if (!std::isfinite(u)) return 0.5;
  if (u < eps) return eps;
  if (u > 1.0 - eps) return 1.0 - eps;

  return u;
}


inline double inv_logit_cpp(double x) {

  if (x >= 0.0) {
    const double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  }

  const double z = std::exp(x);

  return z / (1.0 + z);
}


inline double pnorm_cpp(double x) {

  return 0.5 * std::erfc(-x / std::sqrt(2.0));
}


inline double dnorm_cpp(double x) {

  static const double inv_sqrt_2pi = 0.39894228040143267794;

  return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}


inline double qnorm_cpp(double u) {

  u = clamp01_cpp(u);

  try {

    boost::math::normal_distribution<double> nd(0.0, 1.0);

    return boost::math::quantile(nd, u);

  } catch (...) {

    if (u <= 0.0) return -37.0;
    if (u >= 1.0) return 37.0;

    return 0.0;
  }
}


inline double pweibull_cpp(double y,
                           double shape,
                           double scale) {

  if (y <= 0.0) return 0.0;
  if (shape <= 0.0 || scale <= 0.0) return NAN;

  const double ys = y / scale;

  if (ys <= 0.0 || !std::isfinite(ys)) return NAN;

  const double z = std::pow(ys, shape);

  if (!std::isfinite(z)) {
    return 1.0;
  }

  return -std::expm1(-z);
}


inline double dweibull_cpp(double y,
                           double shape,
                           double scale) {

  if (y <= 0.0) return 0.0;
  if (shape <= 0.0 || scale <= 0.0) return NAN;

  const double ys = y / scale;

  if (ys <= 0.0 || !std::isfinite(ys)) return NAN;

  const double z = std::pow(ys, shape);

  if (!std::isfinite(z)) return 0.0;

  const double logdens =
    std::log(shape) -
    std::log(scale) +
    (shape - 1.0) * std::log(ys) -
    z;

  if (!std::isfinite(logdens)) return 0.0;
  if (logdens < -745.0) return 0.0;
  if (logdens > 700.0) return std::exp(700.0);

  return std::exp(logdens);
}


// ============================================================
// Zero-inflated Weibull marginal functions
// ============================================================

inline double F_ZI_cpp(double y,
                       double pi,
                       double shape,
                       double scale) {

  pi = clamp01_cpp(pi);

  if (y < 0.0) return 1e-10;
  if (y == 0.0) return clamp01_cpp(pi);

  const double H = pweibull_cpp(y, shape, scale);

  if (!std::isfinite(H)) return clamp01_cpp(pi);

  const double out = pi + (1.0 - pi) * H;

  return clamp01_cpp(out);
}


inline double f_ZI_cpp(double y,
                       double pi,
                       double shape,
                       double scale) {

  if (y <= 0.0) return 0.0;

  pi = clamp01_cpp(pi);

  const double dens = dweibull_cpp(y, shape, scale);

  if (!std::isfinite(dens)) return 0.0;

  return (1.0 - pi) * dens;
}


// ============================================================
// Gaussian copula density
// ============================================================

inline double c_gauss_cpp(double u1,
                          double u2,
                          double rho) {

  u1 = clamp01_cpp(u1);
  u2 = clamp01_cpp(u2);

  rho = std::max(-0.999999, std::min(0.999999, rho));

  const double z1 = qnorm_cpp(u1);
  const double z2 = qnorm_cpp(u2);

  const double one_minus_r2 = std::max(1.0 - rho * rho, 1e-12);

  const double log_c =
    -0.5 * std::log(one_minus_r2)
    -0.5 *
    (
        rho * rho * (z1 * z1 + z2 * z2) -
          2.0 * rho * z1 * z2
    ) / one_minus_r2;

  if (!std::isfinite(log_c)) return 0.0;
  if (log_c > 700.0) return std::exp(700.0);
  if (log_c < -745.0) return 0.0;

  return std::exp(log_c);
}


// ============================================================
// Bivariate zero-zero probability
// ============================================================

inline double p00_cpp(double pi1,
                      double pi2,
                      double rho,
                      const double eps) {

  pi1 = clamp01_cpp(pi1);
  pi2 = clamp01_cpp(pi2);

  rho = std::max(-0.999999999, std::min(0.999999999, rho));

  if (std::abs(rho) < 1e-12) {

    double out = pi1 * pi2;

    if (!std::isfinite(out)) return eps;
    if (out < eps) return eps;
    if (out > 1.0) return 1.0;

    return out;
  }

  if (rho > 1.0 - 1e-8) {

    double out = std::min(pi1, pi2);

    if (!std::isfinite(out)) return eps;
    if (out < eps) return eps;
    if (out > 1.0) return 1.0;

    return out;
  }

  if (rho < -1.0 + 1e-8) {

    double out = std::max(pi1 + pi2 - 1.0, 0.0);

    if (!std::isfinite(out)) return eps;
    if (out < eps) return eps;
    if (out > 1.0) return 1.0;

    return out;
  }

  const double one_minus_r2 = std::max(1.0 - rho * rho, 1e-14);
  const double sd = std::sqrt(one_minus_r2);

  double upper_pi;
  double a_cond_pi;

  if (pi2 <= pi1) {
    upper_pi = pi2;
    a_cond_pi = pi1;
  } else {
    upper_pi = pi1;
    a_cond_pi = pi2;
  }

  if (!std::isfinite(upper_pi) || upper_pi <= eps) return eps;

  const double z_upper = qnorm_cpp(upper_pi);
  const double a_cond = qnorm_cpp(a_cond_pi);

  if (!std::isfinite(z_upper) || !std::isfinite(a_cond)) return eps;

  double z_lower = qnorm_cpp(eps);

  if (!std::isfinite(z_lower)) {
    z_lower = -8.0;
  } else {
    z_lower = std::min(z_lower, -8.0);
  }

  if (z_upper <= z_lower) return eps;

  auto integrand = [a_cond, rho, sd](double z) -> double {

    const double arg = (a_cond - rho * z) / sd;
    const double val = pnorm_cpp(arg) * dnorm_cpp(z);

    if (!std::isfinite(val)) return 0.0;

    return val;
  };

  double error = 0.0;
  double L1 = 0.0;

  const int max_depth = 10;
  const double tol = 1e-10;

  double out =
    boost::math::quadrature::gauss_kronrod<double, 31>::integrate(
        integrand,
        z_lower,
        z_upper,
        max_depth,
        tol,
        &error,
        &L1
    );

  if (!std::isfinite(out)) return eps;
  if (out < eps) return eps;
  if (out > 1.0) return 1.0;

  return out;
}


// ============================================================
// Pairwise contributions
// ============================================================

inline double g_0p_cpp(double y2,
                       double pi1,
                       double pi2,
                       double rho,
                       double shape2,
                       double scale2) {

  const double u2 = F_ZI_cpp(y2, pi2, shape2, scale2);

  const double z1 = qnorm_cpp(clamp01_cpp(pi1));
  const double z2 = qnorm_cpp(clamp01_cpp(u2));

  const double sd = std::sqrt(std::max(1.0 - rho * rho, 1e-12));

  const double cond = pnorm_cpp((z1 - rho * z2) / sd);

  return cond * f_ZI_cpp(y2, pi2, shape2, scale2);
}


inline double g_p0_cpp(double y1,
                       double pi1,
                       double pi2,
                       double rho,
                       double shape1,
                       double scale1) {

  const double u1 = F_ZI_cpp(y1, pi1, shape1, scale1);

  const double z1 = qnorm_cpp(clamp01_cpp(u1));
  const double z2 = qnorm_cpp(clamp01_cpp(pi2));

  const double sd = std::sqrt(std::max(1.0 - rho * rho, 1e-12));

  const double cond = pnorm_cpp((z2 - rho * z1) / sd);

  return cond * f_ZI_cpp(y1, pi1, shape1, scale1);
}


inline double g_pp_cpp(double y1,
                       double y2,
                       double pi1,
                       double pi2,
                       double rho,
                       double shape1,
                       double scale1,
                       double shape2,
                       double scale2) {

  const double u1 = F_ZI_cpp(y1, pi1, shape1, scale1);
  const double u2 = F_ZI_cpp(y2, pi2, shape2, scale2);

  const double f1 = f_ZI_cpp(y1, pi1, shape1, scale1);
  const double f2 = f_ZI_cpp(y2, pi2, shape2, scale2);

  return c_gauss_cpp(u1, u2, rho) * f1 * f2;
}


inline double g_pair_cpp(double y1,
                         double y2,
                         double pi1,
                         double pi2,
                         double mu1,
                         double mu2,
                         double d,
                         double shape,
                         double range,
                         double smooth,
                         double gamma_term,
                         double eps) {

  if (range <= 0.0 || shape <= 0.0 || smooth <= 0.0) return eps;

  double rho = rho_matern(d, range, smooth);

  if (!std::isfinite(rho)) return eps;

  rho = std::max(-0.999999999, std::min(0.999999999, rho));

  const double scale1 = mu1 / gamma_term;
  const double scale2 = mu2 / gamma_term;

  if (!std::isfinite(scale1) || !std::isfinite(scale2)) return eps;
  if (scale1 <= 0.0 || scale2 <= 0.0) return eps;

  double val = eps;

  if (y1 == 0.0 && y2 == 0.0) {

    val = p00_cpp(pi1, pi2, rho, eps);

  } else if (y1 == 0.0 && y2 > 0.0) {

    val = g_0p_cpp(y2, pi1, pi2, rho, shape, scale2);

  } else if (y1 > 0.0 && y2 == 0.0) {

    val = g_p0_cpp(y1, pi1, pi2, rho, shape, scale1);

  } else if (y1 > 0.0 && y2 > 0.0) {

    val = g_pp_cpp(
      y1,
      y2,
      pi1,
      pi2,
      rho,
      shape,
      scale1,
      shape,
      scale2
    );
  }

  if (!std::isfinite(val) || val < eps) return eps;

  return val;
}


// ============================================================
// ZI Weibull contribution object
// ============================================================

struct BCLZIWeibullContribution {

  const std::vector<double> beta;
  const std::vector<double> alpha;

  bool has_X;
  bool has_W;

  double shape;
  double range;
  double smooth;
  double gamma_term;
  double eps;

  BCLZIWeibullContribution(const NumericVector& beta_,
                        const NumericVector& alpha_,
                        bool has_X_,
                        bool has_W_,
                        double shape_,
                        double range_,
                        double smooth_,
                        double gamma_term_,
                        double eps_)
    : beta(beta_.begin(), beta_.end()),
      alpha(alpha_.begin(), alpha_.end()),
      has_X(has_X_),
      has_W(has_W_),
      shape(shape_),
      range(range_),
      smooth(smooth_),
      gamma_term(gamma_term_),
      eps(eps_) {

    if (beta.empty()) {
      Rcpp::stop("beta must have positive length.");
    }

    if (alpha.empty()) {
      Rcpp::stop("alpha must have positive length.");
    }

    if (!std::isfinite(shape) || shape <= 0.0) {
      Rcpp::stop("shape must be positive and finite.");
    }

    if (!std::isfinite(range) || range <= 0.0) {
      Rcpp::stop("range must be positive and finite.");
    }

    if (!std::isfinite(smooth) || smooth <= 0.0) {
      Rcpp::stop("smooth must be positive and finite.");
    }

    if (!std::isfinite(gamma_term) || gamma_term <= 0.0) {
      Rcpp::stop("gamma_term must be positive and finite.");
    }
  }

  inline double inv_logit(double x) const {

    return clamp01_cpp(inv_logit_cpp(x));
  }

  inline double operator()(double y1,
                         double y2,
                         double pi1,
                         double pi2,
                         double mu1,
                         double mu2,
                         double d) const {

    const double gij = g_pair_cpp(
      y1,
      y2,
      pi1,
      pi2,
      mu1,
      mu2,
      d,
      shape,
      range,
      smooth,
      gamma_term,
      eps
    );

    return std::log(gij);
  }
};


// ============================================================
// Parameter decoding and evaluation
// ============================================================

inline double bp_zi_eval_ptr_cpp(BCLZIStorage& storage,
                                 const NumericVector& par,
                                 bool has_X,
                                 bool has_W,
                                 double eps,
                                 int nthreads,
                                 int min_pairs_parallel) {

  int idx = 0;

  const double shape = std::exp(par[idx]);
  idx++;

  int p_mu = 1;
  int p_pi = 1;

  if (has_X) {
    p_mu = storage.p_mu[0];
  }

  NumericVector beta(p_mu);

  for (int j = 0; j < p_mu; ++j) {
    beta[j] = par[idx + j];
  }

  idx += p_mu;

  if (has_W) {
    p_pi = storage.p_pi[0];
  }

  NumericVector alpha(p_pi);

  for (int j = 0; j < p_pi; ++j) {
    alpha[j] = par[idx + j];
  }

  idx += p_pi;

  const double range = std::exp(par[idx]);
  idx++;

  const double smooth = std::exp(par[idx]);

  if (!std::isfinite(shape) ||
      !std::isfinite(range) ||
      !std::isfinite(smooth) ||
      shape <= 0.0 ||
      range <= 0.0 ||
      smooth <= 0.0) {
    return 1e20;
  }

  for (int j = 0; j < p_mu; ++j) {
    if (!std::isfinite(beta[j])) return 1e20;
  }

  for (int j = 0; j < p_pi; ++j) {
    if (!std::isfinite(alpha[j])) return 1e20;
  }

  const double gamma_term = std::tgamma(1.0 + 1.0 / shape);

  if (!std::isfinite(gamma_term) || gamma_term <= 0.0) {
    return 1e20;
  }

  BCLZIWeibullContribution contrib(
      beta,
      alpha,
      has_X,
      has_W,
      shape,
      range,
      smooth,
      gamma_term,
      eps
  );

  const double cl = parallel_zi_block_sum(
    storage,
    contrib,
    nthreads,
    min_pairs_parallel
  );

  if (!std::isfinite(cl)) {
    return 1e20;
  }

  return -cl;
}


// ============================================================
// Exported storage constructor
// ============================================================

// [[Rcpp::export]]
SEXP prepare_zi_storage_cpp(List y_list,
                            List pairs_list,
                            List dist_list,
                            Nullable<List> X_list = R_NilValue,
                            Nullable<List> W_list = R_NilValue,
                            bool has_X = false,
                            bool has_W = false) {

  BCLZIStoragePtr ptr(
      new BCLZIStorage(
          y_list,
          pairs_list,
          dist_list,
          X_list,
          W_list,
          has_X,
          has_W
      ),
      true
  );

  return ptr;
}


// ============================================================
// Exported objective using persistent C++ storage
// ============================================================

// [[Rcpp::export]]
double bp_zi_weibull_ptr(SEXP ptr_,
                         NumericVector par,
                         bool has_X = false,
                         bool has_W = false,
                         double eps = 1e-12,
                         int nthreads = 0,
                         int min_pairs_parallel = 20000) {

  BCLZIStoragePtr ptr(ptr_);

  return bp_zi_eval_ptr_cpp(
    *ptr,
    par,
    has_X,
    has_W,
    eps,
    nthreads,
    min_pairs_parallel
  );
}


// ============================================================
// Compatibility wrapper using R lists
// ============================================================

// [[Rcpp::export]]
double bp_zi_weibull_cpp(const NumericVector& par,
                         const List& y_list,
                         const List& pairs_list,
                         const List& dist_list,
                         const Nullable<List>& X_list = R_NilValue,
                         const Nullable<List>& W_list = R_NilValue,
                         const bool has_X = false,
                         const bool has_W = false,
                         const double eps = 1e-12,
                         const int nthreads = 0,
                         const int min_pairs_parallel = 20000) {

  BCLZIStorage storage(
      y_list,
      pairs_list,
      dist_list,
      X_list,
      W_list,
      has_X,
      has_W
  );

  return bp_zi_eval_ptr_cpp(
    storage,
    par,
    has_X,
    has_W,
    eps,
    nthreads,
    min_pairs_parallel
  );
}
