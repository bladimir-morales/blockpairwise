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
  smooth <- 0.5

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
    start_dist = c(
      shape = 1.5
    ),
    start_cov = c(
      range = range,
      smooth = smooth
    ),
    fixed_param = c(smooth = smooth),
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

  expect_named(fit$coefficients_dist, "shape")
  expect_true(all(is.finite(fit$coefficients_dist)))
  expect_true(all(fit$coefficients_dist > 0))

  expect_named(fit$coefficients_cov, c("range", "smooth"))
  expect_true(all(is.finite(fit$coefficients_cov)))
  expect_true(all(fit$coefficients_cov > 0))
  expect_equal(unname(fit$coefficients_cov["smooth"]), smooth, tolerance = 1e-10)

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
  cf_cov <- coef(fit, component = "covariance")
  cf_dist <- coef(fit, component = "distribution")

  expect_true(is.numeric(cf_all))
  expect_true(is.numeric(cf_mu))
  expect_true(is.numeric(cf_zi))
  expect_true(is.numeric(cf_dep))
  expect_true(is.numeric(cf_cov))
  expect_true(is.numeric(cf_dist))

  expect_named(cf_mu, "mu")
  expect_named(cf_zi, c("shape", "mu", "pi"))
  expect_named(cf_dist, "shape")
  expect_named(cf_dep, c("range", "smooth"))
  expect_named(cf_cov, c("range", "smooth"))
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
  smooth <- 0.5

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
    start_dist = c(
      shape = shape
    ),
    start_cov = c(
      range = range,
      smooth = smooth
    ),
    fixed_param = c(smooth = smooth),
    cblocks = "kmeans",
    nblocks = 1,
    fweight = "nn",
    nnk = 3,
    optimizer = "BFGS",
    control = list(maxit = 20)
  )

  expect_s3_class(fit, "bp_zi_fit")
  expect_s3_class(fit, "bp_fit")

  expect_length(fit$coefficients_mu, 2)
  expect_length(fit$coefficients_zi, 2)

  expect_named(fit$coefficients_mu, c("beta0", "beta1"))
  expect_named(fit$coefficients_zi, c("alpha0", "alpha1"))

  expect_true(all(is.finite(fit$coefficients_mu)))
  expect_true(all(is.finite(fit$coefficients_zi)))

  expect_named(fit$coefficients_dist, "shape")
  expect_true(all(is.finite(fit$coefficients_dist)))
  expect_true(all(fit$coefficients_dist > 0))

  expect_named(fit$coefficients_cov, c("range", "smooth"))
  expect_true(all(is.finite(fit$coefficients_cov)))
  expect_true(all(fit$coefficients_cov > 0))
  expect_equal(unname(fit$coefficients_cov["smooth"]), smooth, tolerance = 1e-10)


  cf_mu <- coef(fit, component = "mu")
  cf_pi <- coef(fit, component = "pi")
  cf_zi <- coef(fit, component = "zi")
  cf_dep <- coef(fit, component = "dependence")
  cf_cov <- coef(fit, component = "covariance")
  cf_dist <- coef(fit, component = "distribution")

  expect_named(cf_mu, c("beta0", "beta1"))
  expect_named(cf_pi, c("alpha0", "alpha1"))
  expect_named(cf_zi, c("shape", "beta0", "beta1", "alpha0", "alpha1"))
  expect_named(cf_dist, "shape")
  expect_named(cf_dep, c("range", "smooth"))
  expect_named(cf_cov, c("range", "smooth"))
})


test_that("fit_bp_zi supports exponential case through Matern smooth 0.5", {

  set.seed(123)

  n <- 80

  coords <- cbind(
    runif(n),
    runif(n)
  )

  shape <- 1.5
  mu <- 5
  pi0 <- 0.2
  range <- 0.25
  smooth <- 0.5

  scale <- mu / gamma(1 + 1 / shape)

  z0 <- rbinom(n, 1, pi0) == 1

  y <- numeric(n)
  y[!z0] <- rweibull(
    sum(!z0),
    shape = shape,
    scale = scale
  )

  fit <- fit_bp_zi(
    y = y,
    coords = coords,
    start_beta = log(mu),
    start_alpha = qlogis(pi0),
    start_dist = c(
      shape = shape
    ),
    start_cov = c(
      range = range,
      smooth = smooth
    ),
    fixed_param = c(smooth = smooth),
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

  expect_named(fit$coefficients_dist, "shape")
  expect_named(fit$coefficients_cov, c("range", "smooth"))

  expect_true(all(is.finite(fit$coefficients_dist)))
  expect_true(all(is.finite(fit$coefficients_cov)))

  expect_true(all(fit$coefficients_dist > 0))
  expect_true(all(fit$coefficients_cov > 0))

  expect_equal(unname(fit$coefficients_cov["smooth"]), smooth, tolerance = 1e-10)


  cf_dep <- coef(fit, component = "dependence")
  cf_cov <- coef(fit, component = "covariance")
  cf_zi <- coef(fit, component = "zi")

  expect_named(cf_dep, c("range", "smooth"))
  expect_named(cf_cov, c("range", "smooth"))
  expect_named(cf_zi, c("shape", "mu", "pi"))

  expect_equal(unname(cf_dep["smooth"]), smooth, tolerance = 1e-10)
  expect_equal(unname(cf_cov["smooth"]), smooth, tolerance = 1e-10)
})
