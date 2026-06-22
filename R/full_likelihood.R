full_gaussian_nll <- function(par, data, D = NULL) {

  data <- validate_spatial_data(data)
  validate_par(par)

  n <- nrow(data)

  mu <- rep(par[1], n)
  sigma2 <- par[2]
  phi <- par[3]
  tau2 <- par[4]

  if (is.null(D)) {
    D <- as.matrix(stats::dist(data[, c("x", "y")]))
  }

  Sigma <- tau2 * diag(n) + sigma2 * exp(-D / phi)

  L <- try(chol(Sigma), silent = TRUE)

  if (inherits(L, "try-error")) {
    Sigma <- Sigma + diag(1e-8, n)
    L <- chol(Sigma)
  }

  Xmu <- data$z - mu
  z <- backsolve(L, Xmu, transpose = TRUE)
  logdetSigma <- 2 * sum(log(diag(L)))

  0.5 * (logdetSigma + sum(z^2))
}
