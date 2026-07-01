#' Evaluate Gaussian block-pairwise composite likelihood
#'
#' Evaluates the Gaussian block-pairwise composite log-likelihood
#' for a given regression parameter vector and covariance parameter vector.
#'
#' The mean structure is
#'
#' \deqn{
#'   Y(s_i) = x_i^\top \beta + \varepsilon(s_i),
#' }
#'
#' where the residual spatial process has covariance parameters
#' `sill`, `range`, `nugget`, and `smooth`.
#'
#' @param beta Numeric vector of regression coefficients. Its length must match
#'   the number of columns of the design matrix used in [prepare_bcl_data()].
#'   For an intercept-only model, use a scalar value, e.g. `beta = mean(y)`.
#' @param cov Named numeric vector with covariance parameters:
#'   `sill`, `range`, `nugget`, and `smooth`.
#' @param bcl_data Object created by [prepare_bcl_data()].
#' @param type Character. Either `"marginal"` or `"conditional"`.
#'
#' @return Numeric value of the Gaussian block-pairwise negative composite
#'   log-likelihood.
#'
#' @export
eval_blockpairwise <- function(beta,
                               cov,
                               bcl_data,
                               type = c("marginal", "conditional")) {

  type <- match.arg(type)

  if (missing(beta) || is.null(beta)) {
    stop("beta must be provided.", call. = FALSE)
  }

  if (missing(cov) || is.null(cov)) {
    stop("cov must be provided.", call. = FALSE)
  }

  if (missing(bcl_data) || is.null(bcl_data)) {
    stop("bcl_data must be provided.", call. = FALSE)
  }

  if (!inherits(bcl_data, "bp_precomp")) {
    stop("bcl_data must be an object created by prepare_bcl_data().",
         call. = FALSE)
  }

  beta <- as.numeric(beta)

  if (length(beta) < 1L) {
    stop("beta must contain at least one coefficient.", call. = FALSE)
  }

  if (any(!is.finite(beta))) {
    stop("beta contains non-finite values.", call. = FALSE)
  }

  required_cov <- c("sill", "range", "nugget", "smooth")

  if (is.null(names(cov))) {
    stop(
      "cov must be a named numeric vector with names: ",
      paste(required_cov, collapse = ", "),
      call. = FALSE
    )
  }

  missing_cov <- setdiff(required_cov, names(cov))

  if (length(missing_cov) > 0L) {
    stop(
      "cov is missing the following parameters: ",
      paste(missing_cov, collapse = ", "),
      call. = FALSE
    )
  }

  cov <- as.numeric(cov[required_cov])
  names(cov) <- required_cov

  if (any(!is.finite(cov))) {
    stop("cov contains non-finite values.", call. = FALSE)
  }

  sill   <- unname(cov["sill"])
  range  <- unname(cov["range"])
  nugget <- unname(cov["nugget"])
  smooth <- unname(cov["smooth"])

  if (sill <= 0) {
    stop("cov['sill'] must be positive.", call. = FALSE)
  }

  if (range <= 0) {
    stop("cov['range'] must be positive.", call. = FALSE)
  }

  if (nugget < 0) {
    stop("cov['nugget'] must be non-negative.", call. = FALSE)
  }

  if (smooth <= 0) {
    stop("cov['smooth'] must be positive.", call. = FALSE)
  }

  if (!is.null(bcl_data$ptr)) {
    ptr <- bcl_data$ptr
  } else {
    ptr <- prepare_bcl_ptr(bcl_data)
  }

  t <- sill + nugget
  t2 <- t * t

  if (!is.finite(t) || !is.finite(t2) || t <= 0 || t2 <= 0) {
    stop("Invalid covariance parameters: sill + nugget must be positive.",
         call. = FALSE)
  }

  val <- switch(
    type,

    marginal = bcl_gaussian_marginal_ptr(
      ptr,
      beta = beta,
      sill = sill,
      range = range,
      t = t,
      t2 = t2,
      smooth = smooth
    ),

    conditional = bcl_gaussian_conditional_ptr(
      ptr,
      beta = beta,
      sill = sill,
      range = range,
      t = t,
      t2 = t2,
      smooth = smooth
    )
  )

  if (!is.finite(val)) {
    stop("The block-pairwise likelihood evaluation returned a non-finite value.",
         call. = FALSE)
  }

  n_pairs <- get_bcl_n_pairs(bcl_data)

  dim_pair <- switch(
    type,
    marginal = 2,
    conditional = 1
  )

  const <- n_pairs * dim_pair * log(2 * pi)

  loglik <- -(val + 0.5 * const)

  if (!is.finite(loglik)) {
    stop("The block-pairwise likelihood evaluation returned a non-finite value.",
         call. = FALSE)
  }

  loglik
}


get_bcl_n_pairs <- function(bcl_data) {

  if (!is.null(bcl_data$n_pairs)) {
    return(sum(as.numeric(bcl_data$n_pairs)))
  }

  if (!is.null(bcl_data$valid_pairs)) {

    n_pairs <- vapply(
      bcl_data$valid_pairs,
      function(x) {
        if (is.null(x)) {
          return(0L)
        }

        if (is.matrix(x)) {
          return(nrow(x))
        }

        if (is.data.frame(x)) {
          return(nrow(x))
        }

        if (is.list(x) && all(c("i", "j") %in% names(x))) {
          return(length(x$i))
        }

        if (is.vector(x)) {
          return(length(x) / 2L)
        }

        stop("Unable to infer the number of pairs from valid_pairs.",
             call. = FALSE)
      },
      numeric(1L)
    )

    return(sum(n_pairs))
  }

  if (!is.null(bcl_data$dist_pairs)) {
    n_pairs <- vapply(
      bcl_data$dist_pairs,
      length,
      numeric(1L)
    )

    return(sum(n_pairs))
  }

  stop(
    "Unable to determine the number of block pairs. ",
    "bcl_data must contain either n_pairs, valid_pairs, or dist_pairs.",
    call. = FALSE
  )
}
