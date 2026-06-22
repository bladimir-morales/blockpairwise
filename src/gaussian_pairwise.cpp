// [[Rcpp::depends(BH, RcppParallel)]]
#include <Rcpp.h>
#include <cmath>
#include <limits>

#include "bcl_storage.h"
#include "parallel_bcl.h"
#include "rho_matern.h"

using namespace Rcpp;

struct GaussianMarginalContribution {

  double mu;
  double sigma2;
  double phi;
  double t;
  double t2;
  double nu;

  GaussianMarginalContribution(double mu_,
                               double sigma2_,
                               double phi_,
                               double t_,
                               double t2_,
                               double nu_)
    : mu(mu_), sigma2(sigma2_), phi(phi_),
      t(t_), t2(t2_), nu(nu_) {

    if (t_ <= 0.0 || !std::isfinite(t_)) {
      Rcpp::stop("t must be positive and finite.");
    }

    if (sigma2_ <= 0.0 || !std::isfinite(sigma2_)) {
      Rcpp::stop("sigma2 must be positive and finite.");
    }

    if (phi_ <= 0.0 || !std::isfinite(phi_)) {
      Rcpp::stop("phi must be positive and finite.");
    }

    if (nu_ <= 0.0 || !std::isfinite(nu_)) {
      Rcpp::stop("nu must be positive and finite.");
    }

    if (t2_ <= 0.0 || !std::isfinite(t2_)) {
      Rcpp::stop("t2 must be positive and finite.");
    }
  }

  inline double operator()(double Zi, double Zj, double h) const {

    const double Za = Zi - mu;
    const double Zb = Zj - mu;

    const double rho = rho_matern(h, phi, nu);
    const double r   = rho * sigma2;
    const double s   = t2 - r * r;

    if (s <= 0.0 || !std::isfinite(s)) return 0.0;

    return std::log(s) +
      (t * (Za * Za + Zb * Zb) - 2.0 * r * Za * Zb) / s;
  }
};

struct GaussianConditionalContribution {

  double mu;
  double sigma2;
  double phi;
  double t;
  double t2;
  double nu;
  double logt;

  GaussianConditionalContribution(double mu_,
                                  double sigma2_,
                                  double phi_,
                                  double t_,
                                  double t2_,
                                  double nu_)
    : mu(mu_), sigma2(sigma2_), phi(phi_),
      t(t_), t2(t2_), nu(nu_), logt(0.0) {

    if (t_ <= 0.0 || !std::isfinite(t_)) {
      Rcpp::stop("t must be positive and finite.");
    }

    if (sigma2_ <= 0.0 || !std::isfinite(sigma2_)) {
      Rcpp::stop("sigma2 must be positive and finite.");
    }

    if (phi_ <= 0.0 || !std::isfinite(phi_)) {
      Rcpp::stop("phi must be positive and finite.");
    }

    if (nu_ <= 0.0 || !std::isfinite(nu_)) {
      Rcpp::stop("nu must be positive and finite.");
    }

    if (t2_ <= 0.0 || !std::isfinite(t2_)) {
      Rcpp::stop("t2 must be positive and finite.");
    }

    logt = std::log(t_);
  }

  inline double operator()(double Zi, double Zj, double h) const {

    const double Za = Zi - mu;
    const double Zb = Zj - mu;

    const double rho = rho_matern(h, phi, nu);
    const double r   = rho * sigma2;
    const double s   = t2 - r * r;

    if (s <= 0.0 || !std::isfinite(s)) return 0.0;

    const double diff = Za - (r * Zb / t);

    return std::log(s) - logt + (t / s) * diff * diff;
  }
};


// [[Rcpp::export]]
SEXP prepare_bcl_storage(Rcpp::List Z_list_R,
                         Rcpp::List valid_pairs_list_R,
                         Rcpp::List dist_pairs_list_R) {
  BCLPtr ptr(
      new BCLStorage(Z_list_R, valid_pairs_list_R, dist_pairs_list_R),
      true
  );

  return ptr;
}

// [[Rcpp::export]]
double bcl_gaussian_marginal_ptr(SEXP ptr_,
                                 double mu,
                                 double sigma2,
                                 double phi,
                                 double t,
                                 double t2,
                                 double nu) {

  BCLPtr ptr(ptr_);

  GaussianMarginalContribution contrib(
      mu, sigma2, phi, t, t2, nu
  );

  return 0.5 * parallel_block_pairwise_sum(*ptr, contrib);
}

// [[Rcpp::export]]
double bcl_gaussian_conditional_ptr(SEXP ptr_,
                                    double mu,
                                    double sigma2,
                                    double phi,
                                    double t,
                                    double t2,
                                    double nu) {

  BCLPtr ptr(ptr_);

  GaussianConditionalContribution contrib(
      mu, sigma2, phi, t, t2, nu
  );

  return 0.5 * parallel_block_pairwise_sum(*ptr, contrib);
}
