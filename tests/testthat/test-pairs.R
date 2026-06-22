test_that("assemble_pairs returns consistent lengths", {

  nn_idx <- matrix(c(
    1, 2, 3,
    2, 1, 3,
    3, 1, 2
  ), nrow = 3, byrow = TRUE)

  nn_dists <- matrix(c(
    0, 0.1, 0.2,
    0, 0.1, 0.3,
    0, 0.2, 0.3
  ), nrow = 3, byrow = TRUE)

  out <- assemble_pairs(nn_idx, nn_dists, nnk = 2)

  expect_true(is.list(out))
  expect_equal(ncol(out$valid_pairs), 2)
  expect_equal(nrow(out$valid_pairs), length(out$dist_pairs))
})
