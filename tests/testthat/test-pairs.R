test_that("assemble_pairs returns consistent lengths and no self-pairs", {

  nn_idx <- matrix(c(
    2, 1, 3,
    1, 2, 3,
    1, 3, 2
  ), nrow = 3, byrow = TRUE)

  nn_dists <- matrix(c(
    0.10, 0.00, 0.20,
    0.10, 0.00, 0.30,
    0.20, 0.00, 0.30
  ), nrow = 3, byrow = TRUE)

  out <- assemble_pairs(nn_idx, nn_dists, nnk = 2L)

  expect_true(is.list(out))
  expect_true(is.matrix(out$valid_pairs))
  expect_equal(ncol(out$valid_pairs), 2L)
  expect_equal(nrow(out$valid_pairs), length(out$dist_pairs))

  expect_equal(nrow(out$valid_pairs), 6L)

  expect_false(any(out$valid_pairs[, 1L] == out$valid_pairs[, 2L]))

  expect_equal(
    out$valid_pairs,
    matrix(c(
      1, 2,
      1, 3,
      2, 1,
      2, 3,
      3, 1,
      3, 2
    ), ncol = 2L, byrow = TRUE)
  )

  expect_equal(out$dist_pairs, c(0.10, 0.20, 0.10, 0.30, 0.20, 0.30))
})


test_that("assemble_pairs works when self is not in the first column", {

  nn_idx <- matrix(c(
    2, 3, 1,
    1, 3, 2,
    1, 2, 3
  ), nrow = 3, byrow = TRUE)

  nn_dists <- matrix(c(
    0.10, 0.20, 0.00,
    0.10, 0.30, 0.00,
    0.20, 0.30, 0.00
  ), nrow = 3, byrow = TRUE)

  out <- assemble_pairs(nn_idx, nn_dists, nnk = 2L)

  expect_equal(nrow(out$valid_pairs), 6L)
  expect_false(any(out$valid_pairs[, 1L] == out$valid_pairs[, 2L]))

  expect_equal(
    out$valid_pairs,
    matrix(c(
      1, 2,
      1, 3,
      2, 1,
      2, 3,
      3, 1,
      3, 2
    ), ncol = 2L, byrow = TRUE)
  )
})


test_that("assemble_pairs returns at most n * nnk pairs", {

  nn_idx <- matrix(c(
    1, 1, 2,
    2, 2, 1
  ), nrow = 2L, byrow = TRUE)

  nn_dists <- matrix(c(
    0, 0, 0.5,
    0, 0, 0.5
  ), nrow = 2L, byrow = TRUE)

  out <- assemble_pairs(nn_idx, nn_dists, nnk = 2L)

  expect_true(nrow(out$valid_pairs) <= 2L * 2L)
  expect_false(any(out$valid_pairs[, 1L] == out$valid_pairs[, 2L]))
  expect_equal(length(out$dist_pairs), nrow(out$valid_pairs))
})


test_that("assemble_pairs rejects invalid inputs", {

  nn_idx <- matrix(c(
    1, 2,
    2, 1
  ), nrow = 2L, byrow = TRUE)

  nn_dists <- matrix(c(
    0, 0.1,
    0, 0.1
  ), nrow = 2L, byrow = TRUE)

  expect_error(
    assemble_pairs(nn_idx, nn_dists, nnk = 0L),
    "nnk should be positive"
  )

  expect_error(
    assemble_pairs(nn_idx, nn_dists, nnk = 2L),
    "nn_idx should have at least nnk \\+ 1 columns"
  )
})
