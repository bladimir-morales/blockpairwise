#' Fit Gaussian block-pairwise composite likelihood
#'
#' @param data data.frame or matrix with columns x, y, z.
#' @param par numeric vector c(mu, sigma2, phi, tau2, nu).
#' @param lower lower bounds for optimization.
#' @param upper upper bounds for optimization.
#' @param fixed_par optional full parameter vector. If NULL, uses par.
#' @param free_idx indices of parameters to estimate.
#' @param likelihood "marginal" or "conditional".
#' @param cblocks block construction method.
#' @param nblocks number of spatial blocks.
#' @param fweight pair construction method: "nn" or "dist".
#' @param nnk number of nearest neighbours, used if fweight = "nn".
#' @param w distance cut-off, used if fweight = "dist".
#' @param method optimization method.
#' @param precomp optional precomputed block-pairwise data.
#' @param control optional list passed to optim().
#' @param ... additional arguments passed to optim().
#'
#' @return Object of class "bp_fit".
#' @export
fit_blockpairwise <- function(data,
                              par,
                              lower,
                              upper,
                              fixed_par = NULL,
                              free_idx = NULL,
                              likelihood = c("marginal", "conditional"),
                              cblocks = c("regular", "kmeans", "row", "col"),
                              nblocks = 1L,
                              fweight = c("nn", "dist"),
                              nnk = NULL,
                              w = NULL,
                              method = "L-BFGS-B",
                              precomp = NULL,
                              control = list(),
                              ...) {

  likelihood <- match.arg(likelihood)
  likelihood <- tolower(likelihood)
  cblocks <- match.arg(cblocks)
  fweight <- match.arg(fweight)

  data <- validate_spatial_data(data)
  validate_par(par)

  if (is.null(fixed_par)) fixed_par <- par
  if (is.null(free_idx)) free_idx <- seq_along(fixed_par)

  if (length(lower) != length(fixed_par) || length(upper) != length(fixed_par)) {
    stop("lower and upper must have the same length as par/fixed_par.", call. = FALSE)
  }

  if (is.null(precomp)) {
    precomp <- prepare_bcl_data(
      data = data,
      cblocks = cblocks,
      nblocks = nblocks,
      fweight = fweight,
      nnk = nnk,
      w = w
    )
  } else {
    if (!inherits(precomp, "bp_precomp")) {
      stop("precomp must be an object created by prepare_bcl_data().", call. = FALSE)
    }
  }

  cblocks <- precomp$cblocks
  nblocks <- precomp$nblocks
  fweight <- precomp$fweight
  nnk <- precomp$nnk
  w <- precomp$w

  storage_ptr <- prepare_bcl_ptr(precomp)

  objective <- make_objective(
    fixed_par = fixed_par,
    free_idx = free_idx,
    storage_ptr = storage_ptr,
    likelihood = likelihood
  )

  time_fit <- system.time({
    opt <- stats::optim(
      par = fixed_par[free_idx],
      fn = objective,
      method = method,
      lower = lower[free_idx],
      upper = upper[free_idx],
      hessian = FALSE,
      control = control,
      ...
    )
  })

  par_hat <- fixed_par
  par_hat[free_idx] <- opt$par
  names(par_hat) <- c("mu", "sigma2", "phi", "tau2", "nu")

  out <- list(
    call = match.call(),
    par = par_hat,
    opt = opt,
    value = opt$value,
    convergence = opt$convergence,
    message = opt$message,
    likelihood = likelihood,
    cblocks = cblocks,
    nblocks = nblocks,
    fweight = fweight,
    nnk = nnk,
    w = w,
    time = time_fit[1:3],
    precomp = precomp
  )

  class(out) <- "bp_fit"
  out
}


make_objective <- function(fixed_par, free_idx, storage_ptr, likelihood) {

  force(fixed_par)
  force(free_idx)
  force(storage_ptr)
  force(likelihood)

  function(par_free) {

    par_full <- fixed_par
    par_full[free_idx] <- par_free

    mu <- par_full[1]
    sigma2 <- par_full[2]
    phi <- par_full[3]
    tau2 <- par_full[4]
    nu <- par_full[5]

    if (!is.finite(mu) ||
        !is.finite(sigma2) ||
        !is.finite(phi) ||
        !is.finite(tau2) ||
        !is.finite(nu) ||
        sigma2 <= 0 ||
        phi <= 0 ||
        tau2 < 0 ||
        nu <= 0) {
      return(1e10)
    }

    t <- sigma2 + tau2
    t2 <- t * t

    val <- switch(
      likelihood,
      marginal = bcl_gaussian_marginal_ptr(
        storage_ptr, mu, sigma2, phi, t, t2, nu
      ),
      conditional = bcl_gaussian_conditional_ptr(
        storage_ptr, mu, sigma2, phi, t, t2, nu
      ),
      stop("Unknown likelihood type.", call. = FALSE)
    )

    if (!is.finite(val)) return(1e10)

    val
  }
}
