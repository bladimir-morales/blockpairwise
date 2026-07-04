#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

inline double clamp01_cpp(double u, const double eps = 1e-10) {
  if (!std::isfinite(u)) return 0.5;
  if (u < eps) return eps;
  if (u > 1.0 - eps) return 1.0 - eps;
  return u;
}

inline double inv_logit_cpp(double x) {
  if (x >= 0.0) {
    double z = std::exp(-x);
    return 1.0 / (1.0 + z);
  } else {
    double z = std::exp(x);
    return z / (1.0 + z);
  }
}

inline double F_ZI_cpp(double y,
                       double pi,
                       double shape,
                       double scale) {
  pi = clamp01_cpp(pi);

  if (y < 0.0) return 1e-10;
  if (y == 0.0) return clamp01_cpp(pi);

  double H = R::pweibull(y, shape, scale, 1, 0);
  double out = pi + (1.0 - pi) * H;

  return clamp01_cpp(out);
}

inline double f_ZI_cpp(double y,
                       double pi,
                       double shape,
                       double scale) {
  if (y <= 0.0) return 0.0;

  pi = clamp01_cpp(pi);
  double dens = R::dweibull(y, shape, scale, 0);

  return (1.0 - pi) * dens;
}

inline double c_gauss_cpp(double u1,
                          double u2,
                          double rho) {
  u1 = clamp01_cpp(u1);
  u2 = clamp01_cpp(u2);

  rho = std::max(-0.999999, std::min(0.999999, rho));

  double z1 = R::qnorm(u1, 0.0, 1.0, 1, 0);
  double z2 = R::qnorm(u2, 0.0, 1.0, 1, 0);

  double one_minus_r2 = std::max(1.0 - rho * rho, 1e-12);

  double log_c =
    -0.5 * std::log(one_minus_r2)
    -0.5 * (
        (rho * rho * (z1 * z1 + z2 * z2) - 2.0 * rho * z1 * z2)
    / one_minus_r2
    );

  return std::exp(log_c);
}

inline double p00_cpp(double pi1,
                      double pi2,
                      double rho,
                      Function pmvnorm_fun,
                      const double eps) {

  pi1 = clamp01_cpp(pi1);
  pi2 = clamp01_cpp(pi2);

  rho = std::max(-0.999999, std::min(0.999999, rho));

  double z1 = R::qnorm(pi1, 0.0, 1.0, 1, 0);
  double z2 = R::qnorm(pi2, 0.0, 1.0, 1, 0);

  NumericVector upper = NumericVector::create(z1, z2);
  NumericVector mean = NumericVector::create(0.0, 0.0);

  NumericMatrix sigma(2, 2);
  sigma(0, 0) = 1.0;
  sigma(1, 1) = 1.0;
  sigma(0, 1) = rho;
  sigma(1, 0) = rho;

  NumericVector ans = pmvnorm_fun(
    Named("upper") = upper,
    Named("mean") = mean,
    Named("sigma") = sigma
  );

  double out = ans[0];

  if (!std::isfinite(out)) return eps;
  if (out < eps) return eps;
  if (out > 1.0) return 1.0;

  return out;
}


inline double g_0p_cpp(double y2,
                       double pi1,
                       double pi2,
                       double rho,
                       double shape2,
                       double scale2) {
  double u2 = F_ZI_cpp(y2, pi2, shape2, scale2);

  double z1 = R::qnorm(clamp01_cpp(pi1), 0.0, 1.0, 1, 0);
  double z2 = R::qnorm(clamp01_cpp(u2), 0.0, 1.0, 1, 0);

  double sd = std::sqrt(std::max(1.0 - rho * rho, 1e-12));

  double cond =
    R::pnorm((z1 - rho * z2) / sd, 0.0, 1.0, 1, 0);

  return cond * f_ZI_cpp(y2, pi2, shape2, scale2);
}

inline double g_p0_cpp(double y1,
                       double pi1,
                       double pi2,
                       double rho,
                       double shape1,
                       double scale1) {
  double u1 = F_ZI_cpp(y1, pi1, shape1, scale1);

  double z1 = R::qnorm(clamp01_cpp(u1), 0.0, 1.0, 1, 0);
  double z2 = R::qnorm(clamp01_cpp(pi2), 0.0, 1.0, 1, 0);

  double sd = std::sqrt(std::max(1.0 - rho * rho, 1e-12));

  double cond =
    R::pnorm((z2 - rho * z1) / sd, 0.0, 1.0, 1, 0);

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
  double u1 = F_ZI_cpp(y1, pi1, shape1, scale1);
  double u2 = F_ZI_cpp(y2, pi2, shape2, scale2);

  double f1 = f_ZI_cpp(y1, pi1, shape1, scale1);
  double f2 = f_ZI_cpp(y2, pi2, shape2, scale2);

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
                         double gamma_term,
                         double eps,
                         Function pmvnorm_fun) {
  if (range <= 0.0 || shape <= 0.0) return eps;

  double rho = std::exp(-d / range);
  rho = std::max(-0.999999, std::min(0.999999, rho));

  double scale1 = mu1 / gamma_term;
  double scale2 = mu2 / gamma_term;

  double val = eps;

  if (y1 == 0.0 && y2 == 0.0) {

    val = p00_cpp(pi1, pi2, rho, pmvnorm_fun, eps);

  } else if (y1 == 0.0 && y2 > 0.0) {

    val = g_0p_cpp(y2, pi1, pi2, rho, shape, scale2);

  } else if (y1 > 0.0 && y2 == 0.0) {

    val = g_p0_cpp(y1, pi1, pi2, rho, shape, scale1);

  } else if (y1 > 0.0 && y2 > 0.0) {

    val = g_pp_cpp(
      y1, y2,
      pi1, pi2,
      rho,
      shape, scale1,
      shape, scale2
    );
  }

  if (!std::isfinite(val) || val < eps) return eps;

  return val;
}

// [[Rcpp::export]]
double bp_zi_weibull_cpp(const NumericVector& par,
                         const List& y_list,
                         const List& pairs_list,
                         const List& dist_list,
                         const Nullable<List>& X_list = R_NilValue,
                         const Nullable<List>& W_list = R_NilValue,
                         const bool has_X = false,
                         const bool has_W = false,
                         const double eps = 1e-12) {

  int B = y_list.size();

  if (B == 0) return 1e20;

  int idx = 0;

  double shape = std::exp(par[idx]);
  idx++;

  List Xs;
  List Ws;

  int p_mu = 1;
  int p_pi = 1;

  if (has_X) {
    Xs = List(X_list);
    NumericMatrix X0 = Xs[0];
    p_mu = X0.ncol();
  }

  NumericVector beta(p_mu);
  for (int j = 0; j < p_mu; ++j) {
    beta[j] = par[idx + j];
  }
  idx += p_mu;

  if (has_W) {
    Ws = List(W_list);
    NumericMatrix W0 = Ws[0];
    p_pi = W0.ncol();
  }

  NumericVector alpha(p_pi);
  for (int j = 0; j < p_pi; ++j) {
    alpha[j] = par[idx + j];
  }
  idx += p_pi;

  double range = std::exp(par[idx]);

  if (!std::isfinite(shape) || !std::isfinite(range)) return 1e20;
  if (shape <= 0.0 || range <= 0.0) return 1e20;

  Environment mvtnorm_env = Environment::namespace_env("mvtnorm");
  Function pmvnorm_fun = mvtnorm_env["pmvnorm"];

  double gamma_term = std::tgamma(1.0 + 1.0 / shape);

  if (!std::isfinite(gamma_term) || gamma_term <= 0.0) return 1e20;

  double cl = 0.0;

  for (int b = 0; b < B; ++b) {
    NumericVector yk = y_list[b];
    IntegerMatrix pairs = pairs_list[b];
    NumericVector dist = dist_list[b];

    int n = yk.size();
    int npairs = pairs.nrow();

    if (npairs <= 0) continue;

    NumericVector mu(n);
    NumericVector pi(n);

    if (has_X) {
      NumericMatrix Xk = Xs[b];

      for (int i = 0; i < n; ++i) {
        double eta = 0.0;

        for (int j = 0; j < p_mu; ++j) {
          eta += Xk(i, j) * beta[j];
        }

        mu[i] = std::exp(eta);
      }

    } else {
      double mu0 = std::exp(beta[0]);

      for (int i = 0; i < n; ++i) {
        mu[i] = mu0;
      }
    }

    if (has_W) {
      NumericMatrix Wk = Ws[b];

      for (int i = 0; i < n; ++i) {
        double eta = 0.0;

        for (int j = 0; j < p_pi; ++j) {
          eta += Wk(i, j) * alpha[j];
        }

        pi[i] = clamp01_cpp(inv_logit_cpp(eta));
      }

    } else {
      double pi0 = clamp01_cpp(inv_logit_cpp(alpha[0]));

      for (int i = 0; i < n; ++i) {
        pi[i] = pi0;
      }
    }

    for (int l = 0; l < npairs; ++l) {
      int i = pairs(l, 0) - 1;
      int j = pairs(l, 1) - 1;

      if (i < 0 || j < 0 || i >= n || j >= n) continue;

      double gij = g_pair_cpp(
        yk[i],
          yk[j],
            pi[i],
              pi[j],
                mu[i],
                  mu[j],
                    dist[l],
                        shape,
                        range,
                        gamma_term,
                        eps,
                        pmvnorm_fun
      );

      cl += std::log(gij);
    }
  }

  double out = -cl;

  if (!std::isfinite(out)) return 1e20;

  return out;
}


