#ifndef RHO_MATERN_H
#define RHO_MATERN_H

#include <cmath>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

inline double rho_matern(double d, double phi, double nu) {

  if (phi <= 0.0 || nu <= 0.0)
    return NAN;

  if (d <= 0.0)
    return 1.0;

  double a = d / phi;

  if (a < 1e-12)
    return 1.0;

  // ---- half-integer ν = k + 1/2 ----
  double k_real = nu - 0.5;
  double k_round = std::round(k_real);
  bool is_half_integer = std::abs(k_real - k_round) < 1e-10;
  int k = static_cast<int>(k_round);

  if (is_half_integer){
    if (k == 0)
      return std::exp(-a);
    if (k == 1)
      return std::exp(-a) * (1.0 + a);
    if (k == 2)
      return std::exp(-a) * (1.0 + a + (a * a) / 3.0);
    if (k == 3)
      return std::exp(-a) * (1.0 + a + 0.4 * a * a + (a * a * a) / 15.0);
  }

  // ---- general case (thread-safe) ----

  double log_coef = (1.0 - nu) * std::log(2.0)
    - boost::math::lgamma(nu);

  // Boost Bessel K
  double bess = boost::math::cyl_bessel_k(nu, a);

  if (!std::isfinite(bess) || bess <= 0.0)
    return 0.0;

  double log_besselk = std::log(bess);

  double log_rho = log_coef + nu * std::log(a) + log_besselk;

  return std::exp(log_rho);
}

#endif
