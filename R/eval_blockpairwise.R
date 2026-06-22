#' Evaluate block-pairwise composite likelihood
#'
#' Evaluates the Gaussian block-pairwise composite negative log-likelihood
#' for a given parameter vector.
#'
#' @param par Numeric vector with parameters named `mu`,`sigma2`,`phi`,`tau2`, and `nu`.
#' @param bcl_data Object created by [prepare_bcl_data()].
#' @param type Character. Either `"marginal"` or `"conditional"`.
#'
#' @return Numeric value of the block-pairwise negative log-likelihood.
#'
#' @export
eval_blockpairwise <- function(par,
                               bcl_data,
                               type = c("marginal", "conditional")) {

  type <- match.arg(type)

  if (is.null(bcl_data$ptr)) {
    ptr <- prepare_bcl_ptr(bcl_data)
  } else {
    ptr <- bcl_data$ptr
  }

  required <- c("mu", "sigma2", "phi", "tau2", "nu")
  if (!all(required %in% names(par))) {
    stop(
      "par must be a named numeric vector with names: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }

  validate_par(par)

  mu     <- par[["mu"]]
  sigma2 <- par[["sigma2"]]
  phi    <- par[["phi"]]
  tau2   <- par[["tau2"]]
  nu     <- par[["nu"]]

  t  <- sigma2 + tau2
  t2 <- t * t

  if (type == "marginal") {
    bcl_gaussian_marginal_ptr(
      ptr,
      mu = mu,
      sigma2 = sigma2,
      phi = phi,
      t = t,
      t2 = t2,
      nu = nu
    )
  } else {
    bcl_gaussian_conditional_ptr(
      ptr,
      mu = mu,
      sigma2 = sigma2,
      phi = phi,
      t = t,
      t2 = t2,
      nu = nu
    )
  }
}
