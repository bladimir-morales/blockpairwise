test_that("BCL storage can be created with yk and Xk", {

  precomp <- list(
    yk = list(c(1.0, 2.0, 3.0)),

    Xk = list(
      cbind(
        Intercept = c(1.0, 1.0, 1.0),
        x1 = c(0.2, 0.5, 0.8)
      )
    ),

    valid_pairs = list(
      matrix(c(
        1, 2,
        1, 3
      ), ncol = 2L, byrow = TRUE)
    ),

    dist_pairs = list(c(0.1, 0.2)),

    n_pairs = c(2L)
  )

  class(precomp) <- "bp_precomp"

  ptr <- prepare_bcl_ptr(precomp)

  expect_true(!is.null(ptr))
  expect_equal(typeof(ptr), "externalptr")
})


test_that("BCL storage rejects inconsistent X dimensions", {

  precomp <- list(
    yk = list(c(1.0, 2.0, 3.0)),

    Xk = list(
      matrix(1.0, nrow = 2L, ncol = 1L)
    ),

    valid_pairs = list(
      matrix(c(
        1, 2,
        1, 3
      ), ncol = 2L, byrow = TRUE)
    ),

    dist_pairs = list(c(0.1, 0.2)),

    n_pairs = c(2L)
  )

  class(precomp) <- "bp_precomp"

  expect_error(
    prepare_bcl_ptr(precomp),
    "Each X matrix must have the same number of rows as its y vector"
  )
})


test_that("BCL storage rejects out-of-bounds pair indices", {

  precomp <- list(
    yk = list(c(1.0, 2.0, 3.0)),

    Xk = list(
      matrix(1.0, nrow = 3L, ncol = 1L)
    ),

    valid_pairs = list(
      matrix(c(
        1, 2,
        1, 4
      ), ncol = 2L, byrow = TRUE)
    ),

    dist_pairs = list(c(0.1, 0.2)),

    n_pairs = c(2L)
  )

  class(precomp) <- "bp_precomp"

  expect_error(
    prepare_bcl_ptr(precomp),
    "Pair index out of bounds"
  )
})


test_that("C++ marginal and conditional likelihoods are finite from storage pointer", {

  precomp <- list(
    yk = list(c(1.0, 1.5, 2.0, 2.5)),

    Xk = list(
      cbind(
        Intercept = rep(1.0, 4L),
        x1 = c(0.0, 0.2, 0.4, 0.6)
      )
    ),

    valid_pairs = list(
      matrix(c(
        1, 2,
        1, 3,
        2, 4
      ), ncol = 2L, byrow = TRUE)
    ),

    dist_pairs = list(c(0.1, 0.2, 0.3)),

    n_pairs = c(3L)
  )

  class(precomp) <- "bp_precomp"

  ptr <- prepare_bcl_ptr(precomp)

  beta <- c(1.0, 0.5)
  sill <- 1.0
  range <- 0.4
  nugget <- 0.1
  smooth <- 0.5

  t <- sill + nugget
  t2 <- t * t

  val_marginal <- bcl_gaussian_marginal_ptr(
    ptr,
    beta = beta,
    sill = sill,
    range = range,
    t = t,
    t2 = t2,
    smooth = smooth
  )

  val_conditional <- bcl_gaussian_conditional_ptr(
    ptr,
    beta = beta,
    sill = sill,
    range = range,
    t = t,
    t2 = t2,
    smooth = smooth
  )

  expect_true(is.numeric(val_marginal))
  expect_true(length(val_marginal) == 1L)
  expect_true(is.finite(val_marginal))

  expect_true(is.numeric(val_conditional))
  expect_true(length(val_conditional) == 1L)
  expect_true(is.finite(val_conditional))
})
