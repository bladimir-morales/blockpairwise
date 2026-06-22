#' Creating block labels
#'
#' @param coords location coordinates.
#' @param cblocks block construction method.
#' @param nblocks number of spatial blocks.
#'
#' @return Object of class vector with block labels.
#' @export
make_block <- function(coords,
                       cblocks = c("regular", "kmeans", "row", "col"),
                       nblocks = 1L) {

  cblocks <- match.arg(cblocks)

  nblocks <- as.integer(nblocks)

  if (length(nblocks) != 1L || nblocks < 1L) {
    stop("nblocks must be a positive integer.", call. = FALSE)
  }

  coords <- as.matrix(coords)

  if (!is.numeric(coords)) {
    stop("coords must be numeric.", call. = FALSE)
  }

  if (nrow(coords) < 1L) {
    stop("coords must contain at least one row.", call. = FALSE)
  }

  if (ncol(coords) < 2L) {
    stop("coords must have at least two columns.", call. = FALSE)
  }

  if (any(!is.finite(coords))) {
    stop("coords contains non-finite values.", call. = FALSE)
  }

  n <- nrow(coords)

  if (nblocks == 1L) {
    return(factor(rep(1L, n)))
  }

  if (cblocks == "kmeans") {

    if (n < nblocks) {
      stop("nblocks cannot be larger than the number of observations for kmeans.",
           call. = FALSE)
    }

    if (n < 5e4) {

      category <- stats::kmeans(
        coords,
        centers = nblocks,
        algorithm = "Hartigan-Wong"
      )$cluster

    } else {

      m <- min(n, max(5000L, 20L * nblocks))
      sample_idx <- sample.int(n, size = m)

      init_km <- stats::kmeans(
        coords[sample_idx, , drop = FALSE],
        centers = nblocks,
        nstart = 50,
        iter.max = 200,
        algorithm = "Hartigan-Wong"
      )

      category <- stats::kmeans(
        coords,
        centers = init_km$centers,
        algorithm = "Hartigan-Wong",
        nstart = 1,
        iter.max = 200
      )$cluster
    }

    category <- pmin(pmax(as.integer(category), 1L), nblocks)

    return(factor(category, levels = seq_len(nblocks)))
  }

  if (cblocks == "regular") {

    dims <- bp_regular_grid_dims(nblocks)

    nx <- dims["nx"]
    ny <- dims["ny"]

    cutx <- seq(
      min(coords[, 1]),
      max(coords[, 1]),
      length.out = nx + 1L
    )

    cuty <- seq(
      min(coords[, 2]),
      max(coords[, 2]),
      length.out = ny + 1L
    )

    ix <- bp_find_bin(coords[, 1], cutx, nx)
    iy <- bp_find_bin(coords[, 2], cuty, ny)

    category <- (iy - 1L) * nx + ix

    return(factor(category, levels = seq_len(nblocks)))
  }

  if (cblocks == "row") {

    cuty <- seq(
      min(coords[, 2]),
      max(coords[, 2]),
      length.out = nblocks + 1L
    )

    category <- bp_find_bin(coords[, 2], cuty, nblocks)

    return(factor(category, levels = seq_len(nblocks)))
  }

  if (cblocks == "col") {

    cutx <- seq(
      min(coords[, 1]),
      max(coords[, 1]),
      length.out = nblocks + 1L
    )

    category <- bp_find_bin(coords[, 1], cutx, nblocks)

    return(factor(category, levels = seq_len(nblocks)))
  }
}

bp_regular_grid_dims <- function(nblocks) {

  nblocks <- as.integer(nblocks)

  if (length(nblocks) != 1L || nblocks < 1L) {
    stop("nblocks must be a positive integer.", call. = FALSE)
  }

  divs <- which(nblocks %% seq_len(nblocks) == 0L)

  nx <- divs[which.min(abs(divs - sqrt(nblocks)))]
  ny <- nblocks / nx

  c(nx = as.integer(nx), ny = as.integer(ny))
}

bp_find_bin <- function(x, breaks, nbins) {

  out <- findInterval(
    x,
    breaks,
    rightmost.closed = TRUE,
    all.inside = TRUE
  )

  out <- pmin(pmax(out, 1L), nbins)

  as.integer(out)
}
