validate_spatial_data <- function(data) {

  data <- as.data.frame(data)

  if (ncol(data) < 3L) {
    stop("data must have at least three columns: x, y, z.", call. = FALSE)
  }

  if (!all(c("x", "y", "z") %in% names(data))) {
    names(data)[1:3] <- c("x", "y", "z")
  }

  data <- data[, c("x", "y", "z"), drop = FALSE]

  if (!all(vapply(data, is.numeric, logical(1)))) {
    stop("x, y and z must be numeric.", call. = FALSE)
  }

  if (any(!stats::complete.cases(data))) {
    stop("data contains missing values.", call. = FALSE)
  }

  data
}


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
