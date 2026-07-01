#' Print method for bp_fit objects
#'
#' @param x object of class bp_fit.
#' @param ... ignored.
#' @export
print.bp_fit <- function(x, ...) {
  cat("Block-pairwise Gaussian fit\n")
  cat("Likelihood:", x$likelihood, "\n")
  cat("Blocks:", x$nblocks, "using", x$cblocks, "\n")
  cat("Pair rule:", x$nnk,x$fweight, "\n\n")
  cat("Estimators:", "\n\n")
  print(x$par)
  cat("\nBlock-pairwise log-likelihood:", x$loglik, "\n")
  cat("Convergence code:", x$convergence, "\n")
  invisible(x)
}


#' Summary method for bp_fit objects
#'
#' @param object object of class bp_fit.
#' @param ... ignored.
#' @export
summary.bp_fit <- function(object, ...) {
  out <- list(
    par = object$par,
    loglik = object$loglik,
    value = object$value,
    convergence = object$convergence,
    message = object$message,
    time = object$time,
    likelihood = object$likelihood,
    nblocks = object$nblocks,
    cblocks = object$cblocks,
    fweight = object$fweight
  )
  class(out) <- "summary.bp_fit"
  out
}


#' Print method for summary.bp_fit objects
#'
#' @param x object of class summary.bp_fit.
#' @param ... ignored.
#' @export
print.summary.bp_fit <- function(x, ...) {
  cat("Summary of block-pairwise Gaussian fit\n\n")
  print(x$par)
  cat("\nBlock-pairwise log-likelihood:", x$loglik, "\n")
  cat("Convergence code:", x$convergence, "\n")
  if (!is.null(x$message)) cat("Message:", x$message, "\n")
  cat("\nElapsed time:\n")
  print(x$time)
  invisible(x)
}

#' Extract model coefficients
#'
#' @param object Object of class `"bp_fit"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return Numeric vector of estimated parameters.
#'
#' @export
coef.bp_fit <- function(object, ...) {
  object$par
}

#' Extract log-likelihood
#'
#' @param object Object of class `"bp_fit"`.
#' @param ... Further arguments, currently ignored.
#'
#' @return Object of class `"logLik"`.
#'
#' @export
logLik.bp_fit <- function(object, ...) {

  val <- object$loglik

  attr(val, "df") <- length(object$par)
  attr(val, "nobs") <- object$nobs %||% NA_integer_

  class(val) <- "logLik"
  val
}

#' Simulate from fitted block-pairwise model
#'
#' @param object Object of class `"bp_fit"`.
#' @param nsim Number of simulations.
#' @param seed Optional random seed.
#' @param ... Further arguments.
#'
#' @importFrom stats simulate
#' @export
simulate.bp_fit <- function(object, nsim = 1, seed = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)

  stop("Simulation is not implemented yet.", call. = FALSE)
}

#' Predict from fitted block-pairwise model
#'
#' @param object Object of class `"bp_fit"`.
#' @param newcoords Matrix of prediction coordinates.
#' @param ... Further arguments.
#'
#' @importFrom stats predict
#' @export
predict.bp_fit <- function(object, newcoords, ...) {
  stop("Prediction is not implemented yet.", call. = FALSE)
}


#' Print method for bp_precomp objects
#'
#' @param x object of class bp_precomp.
#' @param ... ignored.
#' @export
print.bp_precomp <- function(x, ...) {
  cat("Block-pairwise precomputed data\n")
  cat("Blocks:", x$nblocks, "using", x$cblocks, "\n")
  cat("Pair rule:", x$fweight, "\n")

  if (!is.null(x$nnk)) cat("Nearest neighbours:", x$nnk, "\n")
  if (!is.null(x$w)) cat("Distance cutoff:", x$w, "\n")

  cat("Total observations:", nrow(x$data), "\n")
  cat("Total pairs:", sum(x$n_pairs), "\n\n")

  cat("Block sizes:\n")
  print(summary(x$block_sizes))

  cat("\nPairs per block:\n")
  print(summary(x$n_pairs))

  invisible(x)
}

#' Summary method for bp_precomp objects
#'
#' @param object object of class bp_precomp.
#' @param ... ignored.
#' @export
summary.bp_precomp <- function(object, ...) {
  out <- list(
    nblocks = object$nblocks,
    block_sizes = summary(object$block_sizes),
    pairs_per_block = summary(object$n_pairs),
    total_pairs = sum(object$n_pairs),
    distance_summary = summary(unlist(object$dist_pairs)),
    cblocks = object$cblocks,
    fweight = object$fweight
  )

  class(out) <- "summary.bp_precomp"
  out
}

#' Print method for summary.bp_precomp objects
#'
#' @param x object of class summary.bp_precomp.
#' @param ... ignored.
#' @export
print.summary.bp_precomp <- function(x, ...) {
  cat("Summary of block-pairwise precomputed data\n\n")

  cat("Blocks:", x$nblocks, "\n")
  cat("Block method:", x$cblocks, "\n")
  cat("Pair rule:", x$fweight, "\n")
  cat("Total pairs:", x$total_pairs, "\n\n")

  cat("Block sizes:\n")
  print(x$block_sizes)

  cat("\nPairs per block:\n")
  print(x$pairs_per_block)

  cat("\nPairwise distance summary:\n")
  print(x$distance_summary)

  invisible(x)
}
