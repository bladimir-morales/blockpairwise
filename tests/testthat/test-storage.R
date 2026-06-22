test_that("BCL storage can be created", {

  Z_list <- list(c(1, 2, 3))

  valid_pairs_list <- list(
    matrix(c(1, 2,
             1, 3), ncol = 2, byrow = TRUE)
  )

  dist_pairs_list <- list(c(0.1, 0.2))

  ptr <- prepare_bcl_storage(
    Z_list,
    valid_pairs_list,
    dist_pairs_list
  )

  expect_true(!is.null(ptr))
})
