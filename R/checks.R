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

bp_build_design_matrix <- function(X, n) {

  if (is.null(X)) {

    X <- matrix(1.0, nrow = n, ncol = 1L)
    colnames(X) <- "Intercept"
    return(X)
  }

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("X must be numeric.", call. = FALSE)
  }

  if (nrow(X) != n) {
    stop("nrow(X) must be equal to length(y).", call. = FALSE)
  }

  if (ncol(X) < 1L) {
    stop("X must have at least one column.", call. = FALSE)
  }

  if (any(!is.finite(X))) {
    stop("X contains non-finite values.", call. = FALSE)
  }

  if (is.null(colnames(X))) {
    colnames(X) <- paste0("x", seq_len(ncol(X)))
  }

  has_intercept <- any(vapply(
    seq_len(ncol(X)),
    function(j) all(abs(X[, j] - 1.0) < sqrt(.Machine$double.eps)),
    logical(1L)
  ))

  if (!has_intercept) {
    X <- cbind(Intercept = 1.0, X)
  }

  storage.mode(X) <- "double"

  X
}
