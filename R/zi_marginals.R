

bp_zi_marginals <- function(y, mu, shape, pi) {

  n <- length(y)

  if (!(length(pi) %in% c(1L, n)) || any(!is.finite(pi)) ||
      any(pi <= 0 | pi >= 1)) {
    stop("pi must be in (0, 1) and have length 1 or n.")
  }
  if (!(length(mu) %in% c(1L, n)) || any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop("mu must be positive, finite and have length 1 or n.")
  }
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0) {
    stop("shape must be positive and finite.")
  }
  if (any(!is.finite(y)) || any(y < 0)) {
    stop("y must be finite and nonnegative.")
  }

  pos <- which(y > 0)

  a <- stats::qnorm(pi)
  b <- rep(NA_real_, n)
  log1m_pi_h <- rep(NA_real_, n)

  if (length(pos) > 0L) {

    mu_pos <- if (length(mu) == 1L) mu else mu[pos]
    pi_pos <- if (length(pi) == 1L) pi else pi[pos]

    scale <- exp(log(mu_pos) - lgamma(1 + 1 / shape))

    log1m_pi_h[pos] <- log1p(-pi_pos) +
      stats::dweibull(y[pos], shape = shape, scale = scale, log = TRUE)

    b[pos] <- stats::qnorm(log1p(-pi_pos) +
                             stats::pweibull(y[pos],
                                             shape = shape, scale = scale,
                                             lower.tail = FALSE, log.p = TRUE),
                           lower.tail = FALSE, log.p = TRUE)
  }

  list(a = a, b = b, log1m_pi_h = log1m_pi_h)
}



bp_zi_marginals <- function(y, mu, shape, eta_pi) {

  n <- length(y)

  if (!(length(eta_pi) %in% c(1L, n)) ||
      any(!is.finite(eta_pi))) {
    stop("eta_pi must be finite and have length 1 or n.")
  }
  if (!(length(mu) %in% c(1L, n)) || any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop("mu must be positive, finite and have length 1 or n.")
  }
  if (length(shape) != 1L || !is.finite(shape) || shape <= 0) {
    stop("shape must be positive and finite.")
  }
  if (any(!is.finite(y)) || any(y < 0)) {
    stop("y must be finite and nonnegative.")
  }

  # Zero-inflation component
  abs_eta <- abs(eta_pi)
  log1p_exp <- log1p(exp(-abs_eta))

  log1m_pi <- -pmax(eta_pi, 0) - log1p_exp
  log_tail <- -abs_eta - log1p_exp

  a <- sign(eta_pi) *
    stats::qnorm(log_tail, lower.tail = FALSE, log.p = TRUE)

  # Positive component
  pos <- which(y > 0)

  b <- rep(NA_real_, n)
  log1m_pi_h <- rep(NA_real_, n)

  if (length(pos) > 0L) {

    mu_pos <- if (length(mu) == 1L) mu else mu[pos]
    log1m_pi_pos <- if (length(eta_pi) == 1L) {
      log1m_pi
    } else {
      log1m_pi[pos]
    }

    scale <- exp(log(mu_pos) - lgamma(1 + 1 / shape))

    log1m_pi_h[pos] <- log1m_pi_pos +
      stats::dweibull(y[pos], shape = shape,
                      scale = scale, log = TRUE)

    b[pos] <- stats::qnorm(
      log1m_pi_pos + stats::pweibull(y[pos], shape = shape,
                                     scale = scale,
                                     lower.tail = FALSE, log.p = TRUE),
      lower.tail = FALSE, log.p = TRUE)
  }

  list(a = a, b = b, log1m_pi_h = log1m_pi_h)
}




