#' Fit zero-inflated Weibull block-pairwise composite likelihood
#'
#' @param y numeric nonnegative response vector.
#' @param coords numeric matrix/data.frame of spatial coordinates.
#' @param X optional design matrix for the positive Weibull mean. If NULL, intercept-only.
#' @param W optional design matrix for the zero-inflation probability. If NULL, intercept-only.
#' @param start_beta optional initial beta vector for log(mu).
#' @param start_alpha optional initial alpha vector for logit(pi).
#' @param start_zi optional named vector with starts: shape, range.
#' @param fixed_param optional named vector of fixed parameters.
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
                      start_zi = NULL,
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
                      control = list(),
                      ...) {

  optimizer <- match.arg(optimizer)
  cblocks   <- match.arg(cblocks)
  fweight   <- match.arg(fweight)

  # 1. Response validation

  y <- as.numeric(y)

  if (length(y) < 2L) {
    stop("y must contain at least two observations.", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop("y contains non-finite values.", call. = FALSE)
  }

  if (any(y < 0)) {
    stop("y must be nonnegative for the zero-inflated Weibull model.",
         call. = FALSE)
  }

  n <- length(y)

  # 2. Coordinates validation

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
  zi_names <- c("shape", "range")

  all_names <- c("shape", beta_names, alpha_names, "range")


  # 4. Starting values in natural/statistical scale ----

  start_param <- bp_zi_make_start_param(
    y = y,
    coords = coords,
    X = X,
    W = W,
    start_beta = start_beta,
    start_alpha = start_alpha,
    start_zi = start_zi,
    fixed_param = fixed_param,
    beta_names = beta_names,
    alpha_names = alpha_names,
    zi_names = zi_names
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
    alpha_names = alpha_names
  )

  # 6. Optimization scale ----
  # shape and range are optimized on log-scale.
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
    stop("All parameters are fixed. There is nothing to optimize.",
         call. = FALSE)
  }

  par0 <- start_opt[free_names_opt]

  # 7. Bounds ----
  # lower and upper must be named in optimization scale:
  # log_shape, beta*, alpha*, log_range.

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
      stop("precomp must be an object created by prepare_bcl_data().",
           call. = FALSE)
    }
  }

  precomp <- bp_zi_attach_designs(
    precomp = precomp,
    X = X,
    W = W
  )

  # 9. Objective function ----

  objective <- make_objective_bcl_zi(
    start_opt = start_opt,
    free_names_opt = free_names_opt,
    beta_names = beta_names,
    alpha_names = alpha_names,
    precomp = precomp,
    eps = eps
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

  names(beta_hat) <- sub("^beta", "", names(beta_hat))
  names(alpha_hat) <- sub("^alpha", "", names(alpha_hat))

  zi_hat <- par_hat[c("shape", "range")]

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
    names(par_user)[names(par_user) == beta_names] <- "mu"
  }

  if (p_alpha == 1L) {
    par_user[alpha_names] <- coefficients_zi
    names(par_user)[names(par_user) == alpha_names] <- "pi"
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

#' @importFrom mvtnorm pmvnorm
NULL

bp_zi_make_start_param <- function(y,
                                   coords,
                                   X,
                                   W,
                                   start_beta = NULL,
                                   start_alpha = NULL,
                                   start_zi = NULL,
                                   fixed_param = NULL,
                                   beta_names,
                                   alpha_names,
                                   zi_names = c("shape", "range")) {

  p_beta <- ncol(X)
  p_alpha <- ncol(W)

  all_names <- c("shape", beta_names, alpha_names, "range")

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

  zi0 <- c(
    shape = 1.5,
    range = 0.25 * domain_diameter
  )

  if (!is.finite(zi0["range"]) || zi0["range"] <= 0) {
    zi0["range"] <- 1
  }

  if (!is.null(start_zi)) {

    if (is.null(names(start_zi))) {
      stop("start_zi must be a named numeric vector.", call. = FALSE)
    }

    start_zi_names <- names(start_zi)
    start_zi <- as.numeric(start_zi)
    names(start_zi) <- start_zi_names

    if (any(!is.finite(start_zi))) {
      stop("start_zi contains non-finite values.", call. = FALSE)
    }

    unknown_zi <- setdiff(names(start_zi), zi_names)

    if (length(unknown_zi) > 0L) {
      stop(
        "Unknown names in start_zi: ",
        paste(unknown_zi, collapse = ", "),
        ". Valid names are: ",
        paste(zi_names, collapse = ", "),
        call. = FALSE
      )
    }

    zi0[names(start_zi)] <- start_zi
  }

  if (length(fixed_names) > 0L) {

    fixed_beta <- intersect(beta_names, fixed_names)
    fixed_alpha <- intersect(alpha_names, fixed_names)
    fixed_zi <- intersect(zi_names, fixed_names)

    if (length(fixed_beta) > 0L) {
      beta0[fixed_beta] <- fixed_param[fixed_beta]
    }

    if (length(fixed_alpha) > 0L) {
      alpha0[fixed_alpha] <- fixed_param[fixed_alpha]
    }

    if (length(fixed_zi) > 0L) {
      zi0[fixed_zi] <- fixed_param[fixed_zi]
    }
  }

  beta0[is.na(beta0)] <- 0
  alpha0[is.na(alpha0)] <- 0

  out <- c(
    shape = unname(zi0["shape"]),
    beta0,
    alpha0,
    range = unname(zi0["range"])
  )

  out
}

bp_zi_check_param <- function(start_param,
                              beta_names,
                              alpha_names) {

  all_names <- c("shape", beta_names, alpha_names, "range")

  missing_names <- setdiff(all_names, names(start_param))

  if (length(missing_names) > 0L) {
    stop(
      "Missing parameters in start_param: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  shape <- unname(start_param["shape"])
  beta <- start_param[beta_names]
  alpha <- start_param[alpha_names]
  range <- unname(start_param["range"])

  if (!is.finite(shape) || shape <= 0) {
    stop("shape must be positive and finite.", call. = FALSE)
  }

  if (any(!is.finite(beta))) {
    stop("beta coefficients contain non-finite values.", call. = FALSE)
  }

  if (any(!is.finite(alpha))) {
    stop("alpha coefficients contain non-finite values.", call. = FALSE)
  }

  if (!is.finite(range) || range <= 0) {
    stop("range must be positive and finite.", call. = FALSE)
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
    log_range = log(unname(start_param["range"]))
  )
}

bp_zi_map_fixed_to_opt <- function(fixed_param_names) {

  out <- fixed_param_names

  out[out == "shape"] <- "log_shape"
  out[out == "range"] <- "log_range"

  out
}


bp_zi_decode_opt_param <- function(par_opt,
                                   beta_names,
                                   alpha_names) {

  shape <- exp(unname(par_opt["log_shape"]))
  beta <- par_opt[beta_names]
  alpha <- par_opt[alpha_names]
  range <- exp(unname(par_opt["log_range"]))

  par <- c(
    shape = shape,
    beta,
    alpha,
    range = range
  )

  list(
    par = par,
    shape = shape,
    beta = beta,
    alpha = alpha,
    range = range
  )
}

bp_zi_make_cpp_param <- function(par_opt,
                                 beta_names,
                                 alpha_names) {

  c(
    log_shape = unname(par_opt["log_shape"]),
    as.numeric(par_opt[beta_names]),
    as.numeric(par_opt[alpha_names]),
    log_range = unname(par_opt["log_range"])
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
                                  precomp,
                                  eps = 1e-12) {

  force(start_opt)
  force(free_names_opt)
  force(beta_names)
  force(alpha_names)
  force(precomp)
  force(eps)

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
    beta <- decoded$beta
    alpha <- decoded$alpha

    if (!is.finite(shape) ||
        !is.finite(range) ||
        shape <= 0 ||
        range <= 0 ||
        any(!is.finite(beta)) ||
        any(!is.finite(alpha))) {
      return(1e10)
    }

    par_cpp <- bp_zi_make_cpp_param(
      par_opt = par_opt,
      beta_names = beta_names,
      alpha_names = alpha_names
    )

    val <- bp_zi_weibull_cpp(
      par = par_cpp,
      y_list = precomp$yk,
      pairs_list = precomp$valid_pairs,
      dist_list = precomp$dist_pairs,
      X_list = precomp$Xk,
      W_list = precomp$Wk,
      has_X = TRUE,
      has_W = TRUE,
      eps = eps
    )

    if (!is.finite(val)) {
      return(1e10)
    }

    val
  }
}




