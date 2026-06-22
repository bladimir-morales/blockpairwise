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

  n <- nrow(coords)

  if (nblocks == 1L) {
    return(factor(rep(1L, n)))
  }

  if (cblocks == "kmeans") {

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
        coords[sample_idx, ],
        centers = nblocks,
        nstart = 50,
        iter.max = 200,
        algorithm = "Hartigan-Wong"
      )

      category <- stats::kmeans(
        coords,
        centers = init_km$centers[, 1:2],
        algorithm = "Hartigan-Wong",
        nstart = 1,
        iter.max = 200
      )$cluster
    }

    return(factor(category))
  }

  if (cblocks == "regular") {

    n_side <- ceiling(sqrt(nblocks))

    cutx <- seq(min(coords[,1]), max(coords[,1]), length.out = n_side + 1L)
    cuty <- seq(min(coords[,2]), max(coords[,2]), length.out = n_side + 1L)

    category <- interaction(
      findInterval(coords[,1], cutx, rightmost.closed = TRUE),
      findInterval(coords[,2], cuty, rightmost.closed = TRUE),
      sep = "-"
    )

    return(factor(category))
  }

  if (cblocks == "row") {

    cuty <- seq(min(coords[,2]), max(coords[,2]), length.out = nblocks + 1L)

    category <- findInterval(coords[,2], cuty, rightmost.closed = TRUE)
    return(factor(category))
  }

  if (cblocks == "col") {

    cutx <- seq(min(coords[,1]), max(coords[,1]), length.out = nblocks + 1L)

    category <- findInterval(coords[,1], cutx, rightmost.closed = TRUE)
    return(factor(category))
  }
}
