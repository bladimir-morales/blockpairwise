validate_par <- function(par) {

  if (!is.numeric(par) || length(par) != 5L) {
    stop("par must be numeric with length 5: c(mu, sigma2, phi, tau2, nu).",
         call. = FALSE)
  }

  if (is.null(names(par))) {
    names(par) <- c("mu", "sigma2", "phi", "tau2", "nu")
  }

  if (!all(is.finite(par))) {
    stop("par must contain only finite values.", call. = FALSE)
  }

  if (par[2] <= 0 || par[3] <= 0 || par[4] < 0 || par[5] <= 0) {
    stop("Parameters must satisfy sigma2 > 0, phi > 0, tau2 >= 0, nu > 0.",
         call. = FALSE)
  }

  invisible(TRUE)
}
