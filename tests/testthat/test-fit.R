test_that("prepare_bcl_data creates a valid intercept-only precomp object", {

  y <- c(1.1, 0.8, 1.4, 1.2)

  coords <- cbind(
    x = c(0.0, 0.2, 0.4, 0.6),
    y = c(0.0, 0.1, 0.0, 0.1)
  )

  precomp <- prepare_bcl_data(
    y = y,
    coords = coords,
    X = NULL,
    cblocks = "regular",
    nblocks = 1L,
    fweight = "nn",
    nnk = 2L
  )

  expect_s3_class(precomp, "bp_precomp")
  expect_true(!is.null(precomp$yk))
  expect_true(!is.null(precomp$Xk))
  expect_true(!is.null(precomp$valid_pairs))
  expect_true(!is.null(precomp$dist_pairs))

  expect_equal(length(precomp$yk), 1L)
  expect_equal(length(precomp$Xk), 1L)

  expect_equal(ncol(precomp$Xk[[1L]]), 1L)
  expect_true(all(precomp$Xk[[1L]][, 1L] == 1.0))

  expect_equal(length(precomp$dist_pairs[[1L]]), nrow(precomp$valid_pairs[[1L]]))
})


test_that("prepare_bcl_data includes intercept when covariates are provided", {

  y <- c(1.1, 0.8, 1.4, 1.2)

  coords <- cbind(
    x = c(0.0, 0.2, 0.4, 0.6),
    y = c(0.0, 0.1, 0.0, 0.1)
  )

  X <- matrix(
    c(0.0, 0.5, 1.0, 1.5),
    ncol = 1L
  )

  colnames(X) <- "z"

  precomp <- prepare_bcl_data(
    y = y,
    coords = coords,
    X = X,
    cblocks = "regular",
    nblocks = 1L,
    fweight = "nn",
    nnk = 2L
  )

  expect_s3_class(precomp, "bp_precomp")

  expect_equal(ncol(precomp$Xk[[1L]]), 2L)
  expect_true(all(precomp$Xk[[1L]][, 1L] == 1.0))
  expect_equal(precomp$Xk[[1L]][, 2L], as.numeric(X[, 1L]))
})


test_that("eval_blockpairwise returns finite values for marginal and conditional likelihoods", {

  y <- c(1.1, 0.8, 1.4, 1.2)

  coords <- cbind(
    x = c(0.0, 0.2, 0.4, 0.6),
    y = c(0.0, 0.1, 0.0, 0.1)
  )

  precomp <- prepare_bcl_data(
    y = y,
    coords = coords,
    X = NULL,
    cblocks = "regular",
    nblocks = 1L,
    fweight = "nn",
    nnk = 2L
  )

  beta <- c(1.0)

  cov <- c(
    sill = 1.0,
    range = 0.4,
    nugget = 0.1,
    smooth = 0.5
  )

  val_marginal <- eval_blockpairwise(
    beta = beta,
    cov = cov,
    bcl_data = precomp,
    type = "marginal"
  )

  val_conditional <- eval_blockpairwise(
    beta = beta,
    cov = cov,
    bcl_data = precomp,
    type = "conditional"
  )

  expect_true(is.numeric(val_marginal))
  expect_true(length(val_marginal) == 1L)
  expect_true(is.finite(val_marginal))

  expect_true(is.numeric(val_conditional))
  expect_true(length(val_conditional) == 1L)
  expect_true(is.finite(val_conditional))
})


test_that("fit_blockpairwise works for intercept-only model", {

  y <- c(1.1, 0.8, 1.4, 1.2, 1.6, 1.3)

  coords <- cbind(
    x = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    y = c(0.0, 0.1, 0.0, 0.1, 0.0, 0.1)
  )

  fit <- fit_blockpairwise(
    y = y,
    coords = coords,
    X = NULL,
    start_beta = c(mean(y)),
    fixed_param = c(
      sill = 1.0,
      range = 0.5,
      nugget = 0.1,
      smooth = 0.5
    ),
    likelihood = "conditional",
    optimizer = "BFGS",
    par_scale = "natural",
    cblocks = "regular",
    nblocks = 1L,
    fweight = "nn",
    nnk = 2L,
    control = list(maxit = 20L)
  )

  expect_s3_class(fit, "bp_fit")

  expect_true(is.numeric(fit$value))
  expect_true(length(fit$value) == 1L)
  expect_true(is.finite(fit$value))

  expect_true(is.numeric(fit$loglik))
  expect_true(length(fit$loglik) == 1L)
  expect_true(is.finite(fit$loglik))

  expect_equal(fit$p_beta, 1L)
  expect_equal(length(fit$coefficients), 1L)
})


test_that("fit_blockpairwise works with intercept plus covariate", {

  y <- c(1.0, 1.2, 1.5, 1.7, 2.0, 2.2)

  coords <- cbind(
    x = c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    y = c(0.0, 0.1, 0.0, 0.1, 0.0, 0.1)
  )

  X <- matrix(
    c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    ncol = 1L
  )

  colnames(X) <- "z"

  fit <- fit_blockpairwise(
    y = y,
    coords = coords,
    X = X,
    start_beta = c(1.0, 1.0),
    fixed_param = c(
      sill = 1.0,
      range = 0.5,
      nugget = 0.1,
      smooth = 0.5
    ),
    likelihood = "conditional",
    optimizer = "BFGS",
    par_scale = "natural",
    cblocks = "regular",
    nblocks = 1L,
    fweight = "nn",
    nnk = 2L,
    control = list(maxit = 20L)
  )

  expect_s3_class(fit, "bp_fit")

  expect_true(is.finite(fit$value))
  expect_true(is.finite(fit$loglik))

  expect_equal(fit$p_beta, 2L)
  expect_equal(length(fit$coefficients), 2L)

  expect_true("precomp" %in% names(fit))
  expect_s3_class(fit$precomp, "bp_precomp")

  expect_equal(ncol(fit$precomp$Xk[[1L]]), 2L)
  expect_true(all(fit$precomp$Xk[[1L]][, 1L] == 1.0))
})
