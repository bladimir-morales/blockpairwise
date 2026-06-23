#' Fit Gaussian block-pairwise composite likelihood
#'
#' @param y numeric response vector.
#' @param coords numeric matrix/data.frame of spatial coordinates.
#' @param X optional design matrix. If NULL, an intercept-only model is used.
#' @param start_beta optional initial beta vector.
#' @param start_cov optional named vector with covariance starts:
#'   sill, range, nugget, smooth.
#' @param fixed_param optional named vector of fixed parameters.
#' @param likelihood "conditional" or "marginal".
#' @param optimizer optimization method passed to stats::optim().
#' @param par_scale "natural" or "log".
#' @param cblocks block construction method.
#' @param nblocks number of blocks.
#' @param fweight pair construction method: "nn" or "dist".
#' @param nnk number of nearest neighbours if fweight = "nn".
#' @param w distance cut-off if fweight = "dist".
#' @param precomp optional object created by prepare_bcl_data().
#' @param lower lower bounds in optimization scale.
#' @param upper upper bounds in optimization scale.
#' @param eps_nugget small constant used only for log-nugget transformation.
#' @param control list passed to stats::optim().
#' @param ... additional arguments passed to stats::optim().
#'
#' @return Object of class "bp_fit".
#' @export
fit_blockpairwise <- function(y,
                              coords,
                              X = NULL,
                              start_beta = NULL,
                              start_cov = NULL,
                              fixed_param = NULL,
                              likelihood = c("conditional", "marginal"),
                              optimizer = c("L-BFGS-B", "BFGS", "Nelder-Mead", "CG"),
                              par_scale = c("natural", "log"),
                              cblocks = c("regular", "kmeans", "row", "col"),
                              nblocks = 1L,
                              fweight = c("nn", "dist"),
                              nnk = NULL,
                              w = NULL,
                              precomp = NULL,
                              lower = NULL,
                              upper = NULL,
                              eps_nugget = 1e-10,
                              control = list(),
                              ...) {

  likelihood <- match.arg(likelihood)
  optimizer  <- match.arg(optimizer)
  par_scale  <- match.arg(par_scale)
  cblocks    <- match.arg(cblocks)
  fweight    <- match.arg(fweight)

  # 1. Response validation

  y <- as.numeric(y)

  if (length(y) < 2L) {
    stop("y must contain at least two observations.", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop("y contains non-finite values.", call. = FALSE)
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


  # 3. Design matrix validation

  if (is.null(X)) {

    X <- matrix(1.0, nrow = n, ncol = 1L)
    colnames(X) <- "Intercept"

  } else {

    X <- as.matrix(X)

    if (!is.numeric(X)) {
      stop("X must be numeric.", call. = FALSE)
    }

    if (nrow(X) != n) {
      stop("nrow(X) must be equal to length(y).", call. = FALSE)
    }

    if (any(!is.finite(X))) {
      stop("X contains non-finite values.", call. = FALSE)
    }

    if (is.null(colnames(X))) {
      colnames(X) <- paste0("x", seq_len(ncol(X)))
      colnames(X)[1L] <- "Intercept"
    }
  }

  p_beta <- ncol(X)

  beta_names <- paste0("beta_", seq_len(p_beta) - 1L)
  cov_names  <- c("sill", "range", "nugget", "smooth")
  all_names  <- c(beta_names, cov_names)

  # 4. Starting values in natural scale ----

  start_param <- bp_make_start_param(
    y = y,
    coords = coords,
    X = X,
    start_beta = start_beta,
    start_cov = start_cov,
    fixed_param = fixed_param,
    beta_names = beta_names,
    cov_names = cov_names
  )

  start_param <- start_param[all_names]

  # 5. Fixed parameters ----

  if (!is.null(fixed_param)) {

    if (is.null(names(fixed_param))) {
      stop("fixed_param must be a named numeric vector.", call. = FALSE)
    }

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

  bp_check_param(
    start_param = start_param,
    beta_names = beta_names,
    cov_names = cov_names
  )

  # 6. Optimization scale ----

  start_opt <- make_start_opt(
    start_param = start_param,
    beta_names = beta_names,
    par_scale = par_scale,
    eps_nugget = eps_nugget
  )

  fixed_param_names_opt <- map_fixed_to_opt(
    fixed_param_names = fixed_param_names,
    par_scale = par_scale
  )

  all_names_opt <- names(start_opt)
  free_names_opt <- setdiff(all_names_opt, fixed_param_names_opt)

  if (length(free_names_opt) == 0L) {
    stop("All parameters are fixed. There is nothing to optimize.", call. = FALSE)
  }

  par0 <- start_opt[free_names_opt]

  # 7. Bounds ----

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

  storage_ptr <- prepare_bcl_ptr(precomp)

  # 9. Objective function

  objective <- make_objective_bcl(
    start_opt = start_opt,
    free_names_opt = free_names_opt,
    beta_names = beta_names,
    storage_ptr = storage_ptr,
    likelihood = likelihood,
    par_scale = par_scale,
    eps_nugget = eps_nugget
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

  decoded <- decode_opt_param(
    par_opt = par_hat_opt,
    beta_names = beta_names,
    par_scale = par_scale,
    eps_nugget = eps_nugget
  )

  par_hat <- decoded$par
  par_hat <- par_hat[all_names]

  beta_hat <- par_hat[beta_names]
  cov_hat  <- par_hat[cov_names]

  names(beta_hat) <- sub("^beta_", "", names(beta_hat))


  # 12. Output ----

  out <- list(
    call = match.call(),

    par = par_hat,
    par_optim_scale = par_hat_opt,
    par_scale = par_scale,

    coefficients = beta_hat,
    cov_par = cov_hat,

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

    likelihood = likelihood,
    optimizer = optimizer,

    n = n,
    p_beta = p_beta,
    X_colnames = colnames(X),

    cblocks = cblocks,
    nblocks = nblocks,
    fweight = fweight,
    nnk = nnk,
    w = w,
    eps_nugget = eps_nugget,

    time = time_fit[1:3],
    opt = opt,
    precomp = precomp
  )

  class(out) <- "bp_fit"

  out
}


make_start_opt <- function(start_param,
                           beta_names,
                           par_scale = c("natural", "log"),
                           eps_nugget = 1e-10) {

  par_scale <- match.arg(par_scale)

  if (par_scale == "natural") {

    out <- c(
      start_param[beta_names],
      sill   = unname(start_param["sill"]),
      range  = unname(start_param["range"]),
      nugget = unname(start_param["nugget"]),
      smooth = unname(start_param["smooth"])
    )

  } else {

    out <- c(
      start_param[beta_names],
      log_sill   = log(unname(start_param["sill"])),
      log_range  = log(unname(start_param["range"])),
      log_nugget = log(unname(start_param["nugget"]) + eps_nugget),
      log_smooth = log(unname(start_param["smooth"]))
    )
  }

  out
}


map_fixed_to_opt <- function(fixed_param_names,
                             par_scale = c("natural", "log")) {

  par_scale <- match.arg(par_scale)

  out <- fixed_param_names

  if (par_scale == "log") {
    out[out == "sill"]   <- "log_sill"
    out[out == "range"]  <- "log_range"
    out[out == "nugget"] <- "log_nugget"
    out[out == "smooth"] <- "log_smooth"
  }

  out
}


decode_opt_param <- function(par_opt,
                             beta_names,
                             par_scale = c("natural", "log"),
                             eps_nugget = 1e-10) {

  par_scale <- match.arg(par_scale)

  beta <- par_opt[beta_names]

  if (par_scale == "natural") {

    sill   <- unname(par_opt["sill"])
    range  <- unname(par_opt["range"])
    nugget <- unname(par_opt["nugget"])
    smooth <- unname(par_opt["smooth"])

  } else {

    sill   <- exp(unname(par_opt["log_sill"]))
    range  <- exp(unname(par_opt["log_range"]))
    nugget <- exp(unname(par_opt["log_nugget"])) - eps_nugget
    smooth <- exp(unname(par_opt["log_smooth"]))
  }

  if (!is.finite(nugget)) {
    nugget <- NA_real_
  } else {
    nugget <- max(nugget, 0)
  }

  par <- c(
    beta,
    sill = sill,
    range = range,
    nugget = nugget,
    smooth = smooth
  )

  list(
    par = par,
    beta = beta,
    sill = sill,
    range = range,
    nugget = nugget,
    smooth = smooth
  )
}


normalize_bounds <- function(lower,
                             upper,
                             free_names_opt) {

  if (is.null(lower)) {

    lower_opt <- rep(-Inf, length(free_names_opt))
    names(lower_opt) <- free_names_opt

  } else {

    lower_names <- names(lower)

    if (is.null(lower_names)) {
      stop("lower must be a named numeric vector.", call. = FALSE)
    }

    lower <- as.numeric(lower)
    names(lower) <- lower_names

    lower_opt <- lower[free_names_opt]
  }

  if (is.null(upper)) {

    upper_opt <- rep(Inf, length(free_names_opt))
    names(upper_opt) <- free_names_opt

  } else {

    upper_names <- names(upper)

    if (is.null(upper_names)) {
      stop("upper must be a named numeric vector.", call. = FALSE)
    }

    upper <- as.numeric(upper)
    names(upper) <- upper_names

    upper_opt <- upper[free_names_opt]
  }

  if (anyNA(lower_opt)) {
    stop(
      "lower must contain bounds for all free optimization parameters: ",
      paste(free_names_opt, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(upper_opt)) {
    stop(
      "upper must contain bounds for all free optimization parameters: ",
      paste(free_names_opt, collapse = ", "),
      call. = FALSE
    )
  }

  if (any(lower_opt >= upper_opt)) {
    stop("Each lower bound must be strictly smaller than its upper bound.",
         call. = FALSE)
  }

  list(
    lower = lower_opt,
    upper = upper_opt
  )
}


bp_make_start_param <- function(y,
                                coords,
                                X,
                                start_beta = NULL,
                                start_cov = NULL,
                                fixed_param = NULL,
                                beta_names,
                                cov_names = c("sill", "range", "nugget", "smooth")) {

  p_beta <- ncol(X)

  all_names <- c(beta_names, cov_names)


  # 1. Fixed parameters


  if (!is.null(fixed_param)) {

    if (is.null(names(fixed_param))) {
      stop("fixed_param must be a named numeric vector.", call. = FALSE)
    }

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

  free_beta_names <- setdiff(beta_names, fixed_names)
  free_cov_names  <- setdiff(cov_names, fixed_names)


  # 2. Beta starts


  beta0 <- rep(NA_real_, p_beta)
  names(beta0) <- beta_names

  if (!is.null(start_beta)) {

    start_beta <- as.numeric(start_beta)

    if (length(start_beta) != p_beta) {
      stop("length(start_beta) must be equal to ncol(X).",
           call. = FALSE)
    }

    if (any(!is.finite(start_beta))) {
      stop("start_beta contains non-finite values.", call. = FALSE)
    }

    beta0[] <- start_beta

  } else {

    if (length(free_beta_names) > 0L) {

      beta_lm <- tryCatch(
        as.numeric(stats::lm.fit(X, y)$coefficients),
        error = function(e) rep(0, p_beta)
      )

      beta_lm[!is.finite(beta_lm)] <- 0
      beta0[] <- beta_lm
    }
  }

  fixed_beta_names <- intersect(beta_names, fixed_names)

  if (length(fixed_beta_names) > 0L) {
    beta0[fixed_beta_names] <- fixed_param[fixed_beta_names]
  }

  beta0[is.na(beta0)] <- 0


  # 3. Residual variance start


  residuals0 <- as.numeric(y - X %*% beta0)

  v0 <- stats::var(residuals0)

  if (!is.finite(v0) || v0 <= 0) {
    v0 <- stats::var(y)
  }

  if (!is.finite(v0) || v0 <= 0) {
    v0 <- 1
  }


  # 4. Cheap spatial range start


  coord_range <- apply(coords, 2L, function(z) {
    max(z) - min(z)
  })

  domain_diameter <- sqrt(sum(coord_range^2))

  if (!is.finite(domain_diameter) || domain_diameter <= 0) {
    domain_diameter <- 1
  }

  range0 <- 0.25 * domain_diameter

  if (!is.finite(range0) || range0 <= 0) {
    range0 <- 1
  }


  # 5. Covariance starts


  cov0 <- c(
    sill   = NA_real_,
    range  = NA_real_,
    nugget = NA_real_,
    smooth = NA_real_
  )

  if ("sill" %in% free_cov_names) {
    cov0["sill"] <- 0.8 * v0
  }

  if ("range" %in% free_cov_names) {
    cov0["range"] <- range0
  }

  if ("nugget" %in% free_cov_names) {
    cov0["nugget"] <- 0.2 * v0
  }

  if ("smooth" %in% free_cov_names) {
    cov0["smooth"] <- 0.5
  }

  if (!is.null(start_cov)) {

    if (is.null(names(start_cov))) {
      stop("start_cov must be a named numeric vector.", call. = FALSE)
    }

    if (any(!is.finite(start_cov))) {
      stop("start_cov contains non-finite values.", call. = FALSE)
    }

    unknown_cov <- setdiff(names(start_cov), cov_names)

    if (length(unknown_cov) > 0L) {
      stop(
        "Unknown names in start_cov: ",
        paste(unknown_cov, collapse = ", "),
        ". Valid covariance names are: ",
        paste(cov_names, collapse = ", "),
        call. = FALSE
      )
    }

    cov0[names(start_cov)] <- start_cov
  }

  fixed_cov_names <- intersect(cov_names, fixed_names)

  if (length(fixed_cov_names) > 0L) {
    cov0[fixed_cov_names] <- fixed_param[fixed_cov_names]
  }

  if (is.na(cov0["sill"])) {
    cov0["sill"] <- 0.8 * v0
  }

  if (is.na(cov0["range"])) {
    cov0["range"] <- range0
  }

  if (is.na(cov0["nugget"])) {
    cov0["nugget"] <- 0.2 * v0
  }

  if (is.na(cov0["smooth"])) {
    cov0["smooth"] <- 0.5
  }

  c(
    beta0,
    cov0[cov_names]
  )
}


bp_check_param <- function(start_param,
                           beta_names,
                           cov_names) {

  all_names <- c(beta_names, cov_names)

  missing_names <- setdiff(all_names, names(start_param))

  if (length(missing_names) > 0L) {
    stop(
      "Missing parameters in start_param: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }

  beta <- start_param[beta_names]

  sill   <- unname(start_param["sill"])
  range  <- unname(start_param["range"])
  nugget <- unname(start_param["nugget"])
  smooth <- unname(start_param["smooth"])

  if (any(!is.finite(beta))) {
    stop("beta contains non-finite values.", call. = FALSE)
  }

  if (!is.finite(sill) || sill <= 0) {
    stop("sill must be positive and finite.", call. = FALSE)
  }

  if (!is.finite(range) || range <= 0) {
    stop("range must be positive and finite.", call. = FALSE)
  }

  if (!is.finite(nugget) || nugget < 0) {
    stop("nugget must be non-negative and finite.", call. = FALSE)
  }

  if (!is.finite(smooth) || smooth <= 0) {
    stop("smooth must be positive and finite.", call. = FALSE)
  }

  invisible(TRUE)
}


make_objective_bcl <- function(start_opt,
                               free_names_opt,
                               beta_names,
                               storage_ptr,
                               likelihood,
                               par_scale = c("natural", "log"),
                               eps_nugget = 1e-10) {

  par_scale <- match.arg(par_scale)

  force(start_opt)
  force(free_names_opt)
  force(beta_names)
  force(storage_ptr)
  force(likelihood)
  force(par_scale)
  force(eps_nugget)

  function(par_free_opt) {

    par_opt <- start_opt
    par_opt[free_names_opt] <- par_free_opt

    decoded <- decode_opt_param(
      par_opt = par_opt,
      beta_names = beta_names,
      par_scale = par_scale,
      eps_nugget = eps_nugget
    )

    beta   <- decoded$beta
    sill   <- decoded$sill
    range  <- decoded$range
    nugget <- decoded$nugget
    smooth <- decoded$smooth

    if (any(!is.finite(beta)) ||
        !is.finite(sill) ||
        !is.finite(range) ||
        !is.finite(nugget) ||
        !is.finite(smooth) ||
        sill <= 0 ||
        range <= 0 ||
        nugget < 0 ||
        smooth <= 0) {
      return(1e10)
    }

    t <- sill + nugget
    t2 <- t * t

    if (!is.finite(t) || !is.finite(t2) || t <= 0 || t2 <= 0) {
      return(1e10)
    }

    val <- switch(
      likelihood,

      marginal = bcl_gaussian_marginal_ptr(
        storage_ptr,
        as.numeric(beta),
        sill,
        range,
        t,
        t2,
        smooth
      ),

      conditional = bcl_gaussian_conditional_ptr(
        storage_ptr,
        as.numeric(beta),
        sill,
        range,
        t,
        t2,
        smooth
      )
    )

    if (!is.finite(val)) {
      return(1e10)
    }

    val
  }
}
