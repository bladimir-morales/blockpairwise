RMSE <- function(x, par) {
  x <- as.matrix(x)
  sqrt(colMeans((sweep(x, 2L, par, "-"))^2))
}


MSE <- function(x, par) {
  x <- as.matrix(x)
  colMeans((sweep(x, 2L, par, "-"))^2)
}


relative_efficiency <- function(sim, par) {
  sim <- as.matrix(sim)
  p <- length(par)
  n <- nrow(sim)

  dif <- sweep(sim, 2L, par, "-")
  inv_godambe_mc <- crossprod(dif) / n

  re_global <- sqrt(det(inv_godambe_mc))^(1 / p)
  re_param <- sqrt(diag(inv_godambe_mc))

  c(re_param, global = re_global)
}


asymptotic_relative_efficiency <- function(sim, par, inv_fisher_mc) {
  sim <- as.matrix(sim)
  p <- length(par)
  n <- nrow(sim)

  dif <- sweep(sim, 2L, par, "-")
  inv_godambe_mc <- crossprod(dif) / n

  are_global <- (sqrt(det(inv_fisher_mc)) / sqrt(det(inv_godambe_mc)))^(1 / p)
  are_param <- sqrt(diag(inv_fisher_mc)) / sqrt(diag(inv_godambe_mc))

  c(are_param, global = are_global)
}
