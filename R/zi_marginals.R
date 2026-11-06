

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





