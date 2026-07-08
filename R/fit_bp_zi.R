#' Fit zero-inflated Weibull block-pairwise composite likelihood
#'
#' @param y numeric nonnegative response vector.
#' @param coords numeric matrix/data.frame of spatial coordinates.
#' @param X optional design matrix for the positive Weibull mean. If NULL, intercept-only.
#' @param W optional design matrix for the zero-inflation probability. If NULL, intercept-only.
#' @param start_beta optional initial beta vector for log(mu).
#' @param start_alpha optional initial alpha vector for logit(pi).
#' @param start_dist optional named vector with positive distributional starts.
#'   For the zero-inflated Weibull model, currently use start_dist = c(shape = ...).
#' @param start_cov optional named vector with covariance starts: range, smooth.
#'   The exponential correlation model is obtained by fixing smooth = 0.5.
#' @param fixed_param optional named vector of fixed parameters. For example,
#'   use fixed_param = c(smooth = 0.5) to fit the exponential correlation model.
#' @param optimizer optimization method passed to stats::optim().
#' @param cblocks block construction method.
#' @param nblocks number of blocks.
#' @param fweight pair construction method: "nn" or "dist".
#' @param nnk number of nearest neighbours if fweight = "nn".
#' @param w distance cut-off if fweight = "dist".
#' @param precomp optional object created by prepare_bcl_data().
#' @param lower lower bounds in optimization scale.
#' @param upper upper bounds in optimization scale.
#' @param eps numerical lower bound for pairwise contributions.
#' @param nthreads optional number of C++ threads used in the block-pairwise
#'   likelihood evaluation. If NULL, the backend chooses automatically.
#' @param control list passed to stats::optim().
#' @param ... additional arguments passed to stats::optim().
#'
#' @return Object of class "bp_zi_fit".
#' @export
fit_bp_zi <- function(y,
                      coords,
                      X = NULL,
                      W = NULL,
                      start_beta = NULL,
                      start_alpha = NULL,
                      start_dist = NULL,
                      start_cov = NULL,
                      fixed_param = NULL,
                      optimizer = c("L-BFGS-B", "BFGS", "Nelder-Mead", "CG"),
                      cblocks = c("regular", "kmeans", "row", "col"),
                      nblocks = 1L,
                      fweight = c("nn", "dist"),
                      nnk = NULL,
                      w = NULL,
                      precomp = NULL,
                      lower = NULL,
                      upper = NULL,
                      eps = 1e-12,
                      nthreads = NULL,
                      control = list(),
                      ...) {

  optimizer <- match.arg(optimizer)
  cblocks <- match.arg(cblocks)
  fweight <- match.arg(fweight)

  if (is.null(nthreads)) {
    nthreads <- 0L
  } else {
    nthreads <- as.integer(nthreads)

    if (length(nthreads) != 1L || is.na(nthreads) || nthreads < 1L) {
      stop("nthreads must be NULL or a positive integer.", call. = FALSE)
    }
  }

  # 1. Response validation ----

  y <- as.numeric(y)

  if (length(y) < 2L) {
    stop("y must contain at least two observations.", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop("y contains non-finite values.", call. = FALSE)
  }

  if (any(y < 0)) {
    stop(
      "y must be nonnegative for the zero-inflated Weibull model.",
      call. = FALSE
    )
  }

  n <- length(y)

  # 2. Coordinates validation ----

  coords <- as.matrix(coords)

  if (!is.numeric(coords)) {
    stop("coords must be numeric.", call. = FALSE)
  }

  if (nrow(coords) != n) {
    stop("nrow(coords) must be equal to length(y).", call. = FALSE)
  }

  if (ncol(coords) < 1L) {
    stop("coords must have at least one coordinate column.", call. = FALSE)
  }

  if (any(!is.finite(coords))) {
    stop("coords contains non-finite values.", call. = FALSE)
  }

  # 3. Design matrices ----

  X <- bp_build_design_matrix(X, n)
  W <- bp_build_design_matrix(W, n)

  p_beta <- ncol(X)
  p_alpha <- ncol(W)

  beta_names <- paste0("beta", seq_len(p_beta) - 1L)
  alpha_names <- paste0("alpha", seq_len(p_alpha) - 1L)

  dist_names <- c("shape")
  cov_names <- c("range", "smooth")

  all_names <- c(dist_names, beta_names, alpha_names, cov_names)

  # 4. Starting values in natural/statistical scale ----

  start_param <- bp_zi_make_start_param(
    y = y,
    coords = coords,
    X = X,
    W = W,
    start_beta = start_beta,
    start_alpha = start_alpha,
    start_dist = start_dist,
    start_cov = start_cov,
    fixed_param = fixed_param,
    beta_names = beta_names,
    alpha_names = alpha_names,
    dist_names = dist_names,
    cov_names = cov_names
  )

  start_param <- start_param[all_names]

  # 5. Fixed parameters ----

  if (!is.null(fixed_param)) {

    if (is.null(names(fixed_param))) {
      stop("fixed_param must be a named numeric vector.", call. = FALSE)
    }

    fixed_param <- as.numeric(fixed_param) |>
      stats::setNames(names(fixed_param))

    if (any(!is.finite(fixed_param))) {
      stop("fixed_param contains non-finite values.", call. = FALSE)
    }

    unknown_fixed <- setdiff(names(fixed_param), all_names)

    if (length(unknown_fixed) > 0L) {
      stop(
        "Unknown names in fixed_param: ",
        paste(unknown_fixed, collapse = ", "),
        ". Valid names are: ",
        paste(all_names, collapse = ", "),
        call. = FALSE
      )
    }

    start_param[names(fixed_param)] <- fixed_param
    fixed_param_names <- names(fixed_param)

  } else {

    fixed_param_names <- character(0L)
  }

  bp_zi_check_param(
    start_param = start_param,
    beta_names = beta_names,
    alpha_names = alpha_names,
    dist_names = dist_names,
    cov_names = cov_names
  )

  # 6. Optimization scale ----
  # shape, range and smooth are optimized on log-scale.
  # beta and alpha remain unconstrained.

  start_opt <- bp_zi_make_start_opt(
    start_param = start_param,
    beta_names = beta_names,
    alpha_names = alpha_names
  )

  fixed_param_names_opt <- bp_zi_map_fixed_to_opt(fixed_param_names)

  all_names_opt <- names(start_opt)
  free_names_opt <- setdiff(all_names_opt, fixed_param_names_opt)

  if (length(free_names_opt) == 0L) {
    stop(
      "All parameters are fixed. There is nothing to optimize.",
      call. = FALSE
    )
  }

  par0 <- start_opt[free_names_opt]

  # 7. Bounds ----
  # lower and upper must be named in optimization scale:
  # log_shape, beta*, alpha*, log_range, log_smooth.

  bounds <- normalize_bounds(
    lower = lower,
    upper = upper,
    free_names_opt = free_names_opt
  )

  lower_opt <- bounds$lower
  upper_opt <- bounds$upper

  # 8. Precomputation ----

  if (is.null(precomp)) {

    precomp <- prepare_bcl_data(
      y = y,
      coords = coords,
      X = X,
      cblocks = cblocks,
      nblocks = nblocks,
      fweight = fweight,
      nnk = nnk,
      w = w
    )

  } else {

    if (!inherits(precomp, "bp_precomp")) {
      stop(
        "precomp must be an object created by prepare_bcl_data().",
        call. = FALSE
      )
    }
  }

  precomp <- bp_zi_attach_designs(
    precomp = precomp,
    X = X,
    W = W
  )

  zi_storage_ptr <- prepare_zi_storage_cpp(
    y_list = precomp$yk,
    pairs_list = precomp$valid_pairs,
    dist_list = precomp$dist_pairs,
    X_list = precomp$Xk,
    W_list = precomp$Wk,
    has_X = TRUE,
    has_W = TRUE
  )

  # 9. Objective function ----

  objective <- make_objective_bcl_zi(
    start_opt = start_opt,
    free_names_opt = free_names_opt,
    beta_names = beta_names,
    alpha_names = alpha_names,
    zi_storage_ptr = zi_storage_ptr,
    eps = eps,
    nthreads = nthreads
  )

  # 10. Optimization ----

  time_fit <- system.time({

    if (optimizer == "L-BFGS-B") {

      opt <- stats::optim(
        par = par0,
        fn = objective,
        method = optimizer,
        lower = lower_opt,
        upper = upper_opt,
        control = control,
        ...
      )

    } else {

      opt <- stats::optim(
        par = par0,
        fn = objective,
        method = optimizer,
        control = control,
        ...
      )
    }
  })

  # 11. Decode fitted parameters ----

  par_hat_opt <- start_opt
  par_hat_opt[free_names_opt] <- opt$par

  decoded <- bp_zi_decode_opt_param(
    par_opt = par_hat_opt,
    beta_names = beta_names,
    alpha_names = alpha_names
  )

  par_hat <- decoded$par
  par_hat <- par_hat[all_names]

  beta_hat <- par_hat[beta_names]
  alpha_hat <- par_hat[alpha_names]

  dist_hat <- par_hat[dist_names]
  cov_hat <- par_hat[cov_names]

  # kept for backward compatibility with existing print/coef/tests
  zi_hat <- par_hat[c(dist_names, cov_names)]

  if (p_beta == 1L) {
    coefficients_mu <- exp(unname(beta_hat[1]))
    names(coefficients_mu) <- "mu"
  } else {
    coefficients_mu <- beta_hat
  }

  if (p_alpha == 1L) {
    coefficients_zi <- stats::plogis(unname(alpha_hat[1]))
    names(coefficients_zi) <- "pi"
  } else {
    coefficients_zi <- alpha_hat
  }

  par_user <- par_hat

  if (p_beta == 1L) {
    par_user[beta_names] <- coefficients_mu
    names(par_user)[match(beta_names, names(par_user))] <- "mu"
  }

  if (p_alpha == 1L) {
    par_user[alpha_names] <- coefficients_zi
    names(par_user)[match(alpha_names, names(par_user))] <- "pi"
  }

  # 12. Composite log-likelihood ----

  n_pairs <- get_bcl_n_pairs(precomp)

  loglik <- -opt$value

  # 13. Output ----

  out <- list(
    call = match.call(),

    par = par_user,
    par_optim_scale = par_hat_opt,

    coefficients_mu = coefficients_mu,
    coefficients_zi = coefficients_zi,

    dist_par = dist_hat,
    cov_par = cov_hat,

    # backward-compatible alias
    zi_par = zi_hat,

    start_param = start_param,
    start_optim_scale = start_opt,

    fixed_param = fixed_param,
    fixed_param_names = fixed_param_names,
    free_param_names = setdiff(all_names, fixed_param_names),
    free_param_names_opt = free_names_opt,

    lower = lower_opt,
    upper = upper_opt,

    value = opt$value,
    convergence = opt$convergence,
    message = opt$message,

    loglik = loglik,

    optimizer = optimizer,
    nthreads = nthreads,

    n = n,
    n_pairs = n_pairs,

    p_beta = p_beta,
    p_alpha = p_alpha,

    X_colnames = colnames(X),
    W_colnames = colnames(W),

    cblocks = cblocks,
    nblocks = nblocks,
    fweight = fweight,
    nnk = nnk,
    w = w,

    eps = eps,

    time = time_fit[1:3],
    opt = opt,
    precomp = precomp
  )

  class(out) <- c("bp_zi_fit", "bp_fit")

  out
}


bp_zi_make_start_param <- function(y,
                                   coords,
                                   X,
                                   W,
                                   start_beta = NULL,
                                   start_alpha = NULL,
                                   start_dist = NULL,
                                   start_cov = NULL,
                                   fixed_param = NULL,
                                   beta_names,
                                   alpha_names,
                                   dist_names = c("shape"),
                                   cov_names = c("range", "smooth")) {

  p_beta <- ncol(X)
  p_alpha <- ncol(W)

  all_names <- c(dist_names, beta_names, alpha_names, cov_names)

  if (!is.null(fixed_param)) {

    if (is.null(names(fixed_param))) {
      stop("fixed_param must be a named numeric vector.", call. = FALSE)
    }

    fixed_param_names <- names(fixed_param)
    fixed_param <- as.numeric(fixed_param)
    names(fixed_param) <- fixed_param_names

    if (any(!is.finite(fixed_param))) {
      stop("fixed_param contains non-finite values.", call. = FALSE)
    }

    unknown_fixed <- setdiff(names(fixed_param), all_names)

    if (length(unknown_fixed) > 0L) {
      stop(
        "Unknown names in fixed_param: ",
        paste(unknown_fixed, collapse = ", "),
        ". Valid names are: ",
        paste(all_names, collapse = ", "),
        call. = FALSE
      )
    }

    fixed_names <- names(fixed_param)

  } else {

    fixed_names <- character(0L)
  }

  beta0 <- rep(NA_real_, p_beta)
  names(beta0) <- beta_names

  if (!is.null(start_beta)) {

    start_beta <- as.numeric(start_beta)

    if (length(start_beta) != p_beta) {
      stop("length(start_beta) must be equal to ncol(X).", call. = FALSE)
    }

    if (any(!is.finite(start_beta))) {
      stop("start_beta contains non-finite values.", call. = FALSE)
    }

    beta0[] <- start_beta

  } else {

    ind_pos <- y > 0

    if (sum(ind_pos) >= p_beta) {

      yy <- log(pmax(y[ind_pos], .Machine$double.eps))

      beta_lm <- tryCatch(
        as.numeric(stats::lm.fit(X[ind_pos, , drop = FALSE], yy)$coefficients),
        error = function(e) rep(0, p_beta)
      )

      beta_lm[!is.finite(beta_lm)] <- 0
      beta0[] <- beta_lm

    } else {

      mu0 <- mean(y[y > 0])

      if (!is.finite(mu0) || mu0 <= 0) {
        mu0 <- mean(y)
      }

      if (!is.finite(mu0) || mu0 <= 0) {
        mu0 <- 1
      }

      beta0[] <- 0
      beta0[1] <- log(mu0)
    }
  }

  alpha0 <- rep(NA_real_, p_alpha)
  names(alpha0) <- alpha_names

  if (!is.null(start_alpha)) {

    start_alpha <- as.numeric(start_alpha)

    if (length(start_alpha) != p_alpha) {
      stop("length(start_alpha) must be equal to ncol(W).", call. = FALSE)
    }

    if (any(!is.finite(start_alpha))) {
      stop("start_alpha contains non-finite values.", call. = FALSE)
    }

    alpha0[] <- start_alpha

  } else {

    z0 <- as.numeric(y == 0)

    if (length(unique(z0)) == 2L && nrow(W) > p_alpha) {

      alpha_glm <- tryCatch(
        as.numeric(stats::glm.fit(
          x = W,
          y = z0,
          family = stats::binomial()
        )$coefficients),
        error = function(e) rep(0, p_alpha)
      )

      alpha_glm[!is.finite(alpha_glm)] <- 0
      alpha0[] <- alpha_glm

    } else {

      pi0 <- mean(y == 0)
      pi0 <- min(max(pi0, 1e-4), 1 - 1e-4)

      alpha0[] <- 0
      alpha0[1] <- stats::qlogis(pi0)
    }
  }

  coord_range <- apply(coords, 2L, function(z) max(z) - min(z))
  domain_diameter <- sqrt(sum(coord_range^2))

  if (!is.finite(domain_diameter) || domain_diameter <= 0) {
    domain_diameter <- 1
  }

  dist0 <- c(
    shape = 1.5
  )

  cov0 <- c(
    range = 0.25 * domain_diameter,
    smooth = 0.5
  )

  if (!is.finite(cov0["range"]) || cov0["range"] <= 0) {
    cov0["range"] <- 1
  }

  if (!is.null(start_dist)) {

    if (is.null(names(start_dist))) {
      stop("start_dist must be a named numeric vector.", call. = FALSE)
    }

    start_dist_names <- names(start_dist)
    start_dist <- as.numeric(start_dist)
    names(start_dist) <- start_dist_names

    if (any(!is.finite(start_dist))) {
      stop("start_dist contains non-finite values.", call. = FALSE)
    }

    unknown_dist <- setdiff(names(start_dist), dist_names)

    if (length(unknown_dist) > 0L) {
      stop(
        "Unknown names in start_dist: ",
        paste(unknown_dist, collapse = ", "),
        ". Valid names are: ",
        paste(dist_names, collapse = ", "),
        call. = FALSE
      )
    }

    dist0[names(start_dist)] <- start_dist
  }

  if (!is.null(start_cov)) {

    if (is.null(names(start_cov))) {
      stop("start_cov must be a named numeric vector.", call. = FALSE)
    }

    start_cov_names <- names(start_cov)
    start_cov <- as.numeric(start_cov)
    names(start_cov) <- start_cov_names

    if (any(!is.finite(start_cov))) {
      stop("start_cov contains non-finite values.", call. = FALSE)
    }

    unknown_cov <- setdiff(names(start_cov), cov_names)

    if (length(unknown_cov) > 0L) {
      stop(
        "Unknown names in start_cov: ",
        paste(unknown_cov, collapse = ", "),
        ". Valid names are: ",
        paste(cov_names, collapse = ", "),
        call. = FALSE
      )
    }

    cov0[names(start_cov)] <- start_cov
  }

  if (length(fixed_names) > 0L) {

    fixed_beta <- intersect(beta_names, fixed_names)
    fixed_alpha <- intersect(alpha_names, fixed_names)
    fixed_dist <- intersect(dist_names, fixed_names)
    fixed_cov <- intersect(cov_names, fixed_names)

    if (length(fixed_beta) > 0L) {
      beta0[fixed_beta] <- fixed_param[fixed_beta]
    }

    if (length(fixed_alpha) > 0L) {
      alpha0[fixed_alpha] <- fixed_param[fixed_alpha]
    }

    if (length(fixed_dist) > 0L) {
      dist0[fixed_dist] <- fixed_param[fixed_dist]
    }

    if (length(fixed_cov) > 0L) {
      cov0[fixed_cov] <- fixed_param[fixed_cov]
    }
  }

  beta0[is.na(beta0)] <- 0
  alpha0[is.na(alpha0)] <- 0

  out <- c(
    shape = unname(dist0["shape"]),
    beta0,
    alpha0,
    range = unname(cov0["range"]),
    smooth = unname(cov0["smooth"])
  )

  out
}


bp_zi_check_param <- function(start_param,
                              beta_names,
                              alpha_names,
                              dist_names = c("shape"),
                              cov_names = c("range", "smooth")) {

  all_names <- c(dist_names, beta_names, alpha_names, cov_names)

  missing_names <- setdiff(all_names, names(start_param))

  if (length(missing_names) > 0L) {
    stop(
      "Missing parameters in start_param: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  dist <- start_param[dist_names]
  beta <- start_param[beta_names]
  alpha <- start_param[alpha_names]
  cov <- start_param[cov_names]

  if (any(!is.finite(dist)) || any(dist <= 0)) {
    stop(
      "Distribution parameters must be positive and finite.",
      call. = FALSE
    )
  }

  if (any(!is.finite(beta))) {
    stop("beta coefficients contain non-finite values.", call. = FALSE)
  }

  if (any(!is.finite(alpha))) {
    stop("alpha coefficients contain non-finite values.", call. = FALSE)
  }

  if (any(!is.finite(cov)) || any(cov <= 0)) {
    stop(
      "Covariance parameters must be positive and finite.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


bp_zi_make_start_opt <- function(start_param,
                                 beta_names,
                                 alpha_names) {

  c(
    log_shape = log(unname(start_param["shape"])),
    start_param[beta_names],
    start_param[alpha_names],
    log_range = log(unname(start_param["range"])),
    log_smooth = log(unname(start_param["smooth"]))
  )
}


bp_zi_map_fixed_to_opt <- function(fixed_param_names) {

  out <- fixed_param_names

  out[out == "shape"] <- "log_shape"
  out[out == "range"] <- "log_range"
  out[out == "smooth"] <- "log_smooth"

  out
}


bp_zi_decode_opt_param <- function(par_opt,
                                   beta_names,
                                   alpha_names) {

  shape <- exp(unname(par_opt["log_shape"]))
  beta <- par_opt[beta_names]
  alpha <- par_opt[alpha_names]
  range <- exp(unname(par_opt["log_range"]))
  smooth <- exp(unname(par_opt["log_smooth"]))

  par <- c(
    shape = shape,
    beta,
    alpha,
    range = range,
    smooth = smooth
  )

  list(
    par = par,
    shape = shape,
    beta = beta,
    alpha = alpha,
    range = range,
    smooth = smooth
  )
}


bp_zi_make_cpp_param <- function(par_opt,
                                 beta_names,
                                 alpha_names) {

  c(
    log_shape = unname(par_opt["log_shape"]),
    as.numeric(par_opt[beta_names]),
    as.numeric(par_opt[alpha_names]),
    log_range = unname(par_opt["log_range"]),
    log_smooth = unname(par_opt["log_smooth"])
  )
}


bp_zi_attach_designs <- function(precomp,
                                 X,
                                 W) {

  if (is.null(precomp$groups)) {
    stop("precomp must contain groups.", call. = FALSE)
  }

  groups <- precomp$groups

  if (!is.factor(groups)) {
    groups <- factor(groups)
  }

  block_ids <- levels(groups)

  Xk <- vector("list", length(block_ids))
  Wk <- vector("list", length(block_ids))

  names(Xk) <- block_ids
  names(Wk) <- block_ids

  for (b in block_ids) {
    ind <- which(groups == b)

    Xk[[b]] <- X[ind, , drop = FALSE]
    Wk[[b]] <- W[ind, , drop = FALSE]
  }

  precomp$Xk <- Xk
  precomp$Wk <- Wk

  precomp
}


make_objective_bcl_zi <- function(start_opt,
                                  free_names_opt,
                                  beta_names,
                                  alpha_names,
                                  zi_storage_ptr,
                                  eps = 1e-12,
                                  nthreads = 0L) {

  force(start_opt)
  force(free_names_opt)
  force(beta_names)
  force(alpha_names)
  force(zi_storage_ptr)
  force(eps)
  force(nthreads)

  function(par_free_opt) {

    par_opt <- start_opt
    par_opt[free_names_opt] <- par_free_opt

    decoded <- bp_zi_decode_opt_param(
      par_opt = par_opt,
      beta_names = beta_names,
      alpha_names = alpha_names
    )

    shape <- decoded$shape
    range <- decoded$range
    smooth <- decoded$smooth
    beta <- decoded$beta
    alpha <- decoded$alpha

    if (!is.finite(shape) ||
        !is.finite(range) ||
        !is.finite(smooth) ||
        shape <= 0 ||
        range <= 0 ||
        smooth <= 0 ||
        any(!is.finite(beta)) ||
        any(!is.finite(alpha))) {
      return(1e10)
    }

    par_cpp <- bp_zi_make_cpp_param(
      par_opt = par_opt,
      beta_names = beta_names,
      alpha_names = alpha_names
    )

    val <- bp_zi_weibull_ptr(
      ptr_ = zi_storage_ptr,
      par = par_cpp,
      has_X = TRUE,
      has_W = TRUE,
      eps = eps,
      nthreads = nthreads
    )

    if (!is.finite(val)) {
      return(1e10)
    }

    val
  }
}
