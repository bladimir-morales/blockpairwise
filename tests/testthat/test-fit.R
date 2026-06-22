test_that("Gaussian marginal likelihood is finite", {

  Z_list <- list(c(1.1, 0.8, 1.4))

  valid_pairs_list <- list(
    matrix(c(1, 2,
             1, 3,
             2, 3), ncol = 2, byrow = TRUE)
  )

  dist_pairs_list <- list(c(0.1, 0.2, 0.15))

  ptr <- prepare_bcl_storage(
    Z_list,
    valid_pairs_list,
    dist_pairs_list
  )

  val <- bcl_gaussian_marginal_ptr(
    ptr,
    mu = 1,
    sigma2 = 1,
    phi = 0.3,
    t = 1.1,
    t2 = 1.21,
    nu = 0.5
  )

  expect_true(is.finite(val))
})
