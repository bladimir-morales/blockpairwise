#' Print method for zero-inflated block-pairwise fit
#'
#' @param x object of class "bp_zi_fit".
#' @param ... further arguments.
#'
#' @export
print.bp_zi_fit <- function(x, ...) {

  cat("\nZero-inflated Weibull block-pairwise composite likelihood fit\n")
  cat("----------------------------------------------------------------\n")

  cat("Number of observations:", x$n, "\n")
  cat("Number of pairs:", x$n_pairs, "\n")
  cat("Optimizer:", x$optimizer, "\n")
  cat("Convergence:", x$convergence, "\n")

  if (!is.null(x$message)) {
    cat("Message:", x$message, "\n")
  }

  cat("\nComposite log-likelihood:", format(x$loglik, digits = 6), "\n")
  cat("Objective value:", format(x$value, digits = 6), "\n")

  cat("\nZero-inflated marginal parameters:\n")
  print(x$zi_par)

  cat("\nMean component:\n")
  print(x$coefficients_mu)

  cat("\nZero-inflation component:\n")
  print(x$coefficients_zi)

  invisible(x)
}

#' Extract coefficients from zero-inflated block-pairwise fit
#'
#' @param object object of class "bp_zi_fit".
#' @param component one of "all", "mu", "zi", or "dependence".
#' @param ... further arguments.
#'
#' @export
coef.bp_zi_fit <- function(object,
                           component = c("all", "mu", "zi", "dependence"),
                           ...) {

  component <- match.arg(component)

  out <- switch(
    component,

    all = c(
      object$zi_par,
      mu = object$coefficients_mu,
      zi = object$coefficients_zi
    ),

    mu = object$coefficients_mu,

    zi = object$coefficients_zi,

    dependence = object$zi_par
  )

  out
}

#' Log-likelihood method for zero-inflated block-pairwise fit
#'
#' @param object object of class "bp_zi_fit".
#' @param ... further arguments.
#'
#' @export
logLik.bp_zi_fit <- function(object, ...) {

  val <- object$loglik

  attr(val, "df") <- length(object$free_param_names)
  attr(val, "nobs") <- object$n
  class(val) <- "logLik"

  val
}
