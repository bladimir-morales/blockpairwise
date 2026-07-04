test_that("fit_bp_zi returns a valid bp_zi_fit object without covariates", {

  set.seed(123)

  n <- 80

  coords <- cbind(
    runif(n),
    runif(n)
  )

  shape <- 1.6
  mu <- 5
  pi0 <- 0.25
  range <- 0.2

  scale <- mu / gamma(1 + 1 / shape)

  is_zero <- rbinom(n, size = 1, prob = pi0) == 1

  y <- numeric(n)
  y[!is_zero] <- rweibull(
    sum(!is_zero),
    shape = shape,
    scale = scale
  )

  fit <- fit_bp_zi(
    y = y,
    coords = coords,
    start_beta = log(mean(y[y > 0])),
    start_alpha = qlogis(mean(y == 0)),
    start_zi = c(shape = 1.5, range = range),
    cblocks = "kmeans",
    nblocks = 1,
    fweight = "nn",
    nnk = 3,
    optimizer = "BFGS",
    control = list(maxit = 20)
  )

  expect_s3_class(fit, "bp_zi_fit")
  expect_s3_class(fit, "bp_fit")

  expect_true(is.numeric(fit$value))
  expect_true(is.finite(fit$value))

  expect_true(is.numeric(fit$loglik))
  expect_true(is.finite(fit$loglik))

  expect_true(is.numeric(fit$n_pairs))
  expect_gt(fit$n_pairs, 0)

  expect_named(fit$zi_par, c("shape", "range"))
  expect_true(all(is.finite(fit$zi_par)))
  expect_true(all(fit$zi_par > 0))

  expect_named(fit$coefficients_mu, "mu")
  expect_named(fit$coefficients_zi, "pi")

  expect_true(is.finite(fit$coefficients_mu))
  expect_true(fit$coefficients_mu > 0)

  expect_true(is.finite(fit$coefficients_zi))
  expect_true(fit$coefficients_zi > 0)
  expect_true(fit$coefficients_zi < 1)

  expect_true(inherits(logLik(fit), "logLik"))

  cf_all <- coef(fit)
  cf_mu <- coef(fit, component = "mu")
  cf_zi <- coef(fit, component = "zi")
  cf_dep <- coef(fit, component = "dependence")

  expect_true(is.numeric(cf_all))
  expect_true(is.numeric(cf_mu))
  expect_true(is.numeric(cf_zi))
  expect_true(is.numeric(cf_dep))

  expect_named(cf_dep, c("shape", "range"))
})

test_that("fit_bp_zi handles covariates in X and W", {

  set.seed(456)

  n <- 100

  coords <- cbind(
    runif(n),
    runif(n)
  )

  x1 <- rnorm(n)
  w1 <- rnorm(n)

  X <- cbind(1, x1)
  W <- cbind(1, w1)

  colnames(X) <- c("intercept", "x1")
  colnames(W) <- c("intercept", "w1")

  beta <- c(log(4), 0.3)
  alpha <- c(qlogis(0.2), -0.4)

  mu_i <- exp(as.vector(X %*% beta))
  pi_i <- plogis(as.vector(W %*% alpha))

  shape <- 1.5
  range <- 0.2

  scale_i <- mu_i / gamma(1 + 1 / shape)

  is_zero <- rbinom(n, size = 1, prob = pi_i) == 1

  y <- numeric(n)

  y[!is_zero] <- stats::rweibull(
    sum(!is_zero),
    shape = shape,
    scale = scale_i[!is_zero]
  )

  fit <- fit_bp_zi(
    y = y,
    coords = coords,
    X = X,
    W = W,
    start_beta = beta,
    start_alpha = alpha,
    start_zi = c(shape = shape, range = range),
    cblocks = "kmeans",
    nblocks = 1,
    fweight = "nn",
    nnk = 3,
    optimizer = "BFGS",
    control = list(maxit = 20)
  )

  expect_s3_class(fit, "bp_zi_fit")

  expect_length(fit$coefficients_mu, 2)
  expect_length(fit$coefficients_zi, 2)

  expect_named(fit$coefficients_mu, c("0", "1"))
  expect_named(fit$coefficients_zi, c("0", "1"))

  expect_true(all(is.finite(fit$coefficients_mu)))
  expect_true(all(is.finite(fit$coefficients_zi)))

  expect_named(fit$zi_par, c("shape", "range"))
  expect_true(all(fit$zi_par > 0))
})
