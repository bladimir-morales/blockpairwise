#' Print method for zero-inflated block-pairwise fit
#'
#' @param x object of class "bp_zi_fit".
#' @param ... further arguments.
#'
#' @export
print.bp_zi_fit <- function(x, ...) {

  cat("\nZero-inflated Weibull block-pairwise likelihood fit\n")
  cat("-----------------------------------------------------\n")

  cat("Number of observations:", x$n, "\n")
  cat("Number of blocks:", x$nblocks, "using", x$cblocks, "\n")
  cat("Number of pairs:", x$n_pairs, "\n")
  cat("Pair rule:", x$nnk,x$fweight, "\n\n")
  cat("Optimizer:", x$optimizer, "\n")
  cat("Convergence:", x$convergence, "\n")

  cat("\nBlock-pairwise log-likelihood:", format(x$loglik, digits = 6), "\n")

  cat("\nZero-inflation component:\n")
  print(c(x$coefficients_dist, x$coefficients_mu, x$coefficients_zi))

  cat("\nCovariance model:\n")
  print(x$coefficients_cov)

  if (!is.null(x$fixed_param) || !length(x$fixed_param) == 0L) {
    cat("\nFixed parameters:\n")
    print(x$fixed_param)
  }

  cat("\nElapsed time:\n")
  print(x$time)

  invisible(x)
}


#' Extract coefficients from a zero-inflated block-pairwise fit
#'
#' @param object object of class "bp_zi_fit".
#' @param component component to extract.
#' @param ... currently unused.
#'
#' @export
coef.bp_zi_fit <- function(object,
                           component = c(
                             "all",
                             "zero_inflated",
                             "zi",
                             "distribution",
                             "mu",
                             "pi",
                             "covariance",
                             "dependence",
                             "fixed"
                           ),
                           ...) {

  component <- match.arg(component)

  switch(
    component,

    all = c(
      object$coefficients_dist,
      object$coefficients_mu,
      object$coefficients_zi,
      object$coefficients_cov
    ),

    zero_inflated = c(
      object$coefficients_dist,
      object$coefficients_mu,
      object$coefficients_zi
    ),

    zi = c(
      object$coefficients_dist,
      object$coefficients_mu,
      object$coefficients_zi
    ),

    distribution = object$coefficients_dist,

    mu = object$coefficients_mu,

    pi = object$coefficients_zi,

    covariance = object$coefficients_cov,

    dependence = object$coefficients_cov,

    fixed = {
      if (is.null(object$fixed_param) || length(object$fixed_param) == 0L) {
        numeric(0L)
      } else {
        object$fixed_param
      }
    }
  )
}


#' @export
logLik.bp_zi_fit <- function(object, ...) {

  val <- object$loglik

  attr(val, "df") <- length(object$free_param_names)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"

  val
}
