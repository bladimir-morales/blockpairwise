#' Creating block labels
#'
#' @param data data.frame or matrix with columns x, y, z.
#' @param cblocks block construction method.
#' @param nblocks number of spatial blocks.
#'
#' @return Object of class vector with block labels.
#' @export
make_block <- function(data,
                       cblocks = c("regular", "kmeans", "row", "col"),
                       nblocks = 1L) {

  data <- validate_spatial_data(data)
  cblocks <- match.arg(cblocks)

  nblocks <- as.integer(nblocks)
  if (length(nblocks) != 1L || nblocks < 1L) {
    stop("nblocks must be a positive integer.", call. = FALSE)
  }

  n <- nrow(data)

  if (nblocks == 1L) {
    return(factor(rep(1L, n)))
  }

  if (cblocks == "kmeans") {

    if (n < 5e4) {
      category <- stats::kmeans(
        data[, c("x", "y")],
        centers = nblocks,
        algorithm = "Hartigan-Wong"
      )$cluster
    } else {
      m <- min(n, max(5000L, 20L * nblocks))
      sample_idx <- sample.int(n, size = m)

      init_km <- stats::kmeans(
        data[sample_idx, c("x", "y")],
        centers = nblocks,
        nstart = 50,
        iter.max = 200,
        algorithm = "Hartigan-Wong"
      )

      category <- stats::kmeans(
        data[, c("x", "y")],
        centers = init_km$centers[, c("x", "y")],
        algorithm = "Hartigan-Wong",
        nstart = 1,
        iter.max = 200
      )$cluster
    }

    return(factor(category))
  }

  if (cblocks == "regular") {

    n_side <- ceiling(sqrt(nblocks))

    cutx <- seq(min(data$x), max(data$x), length.out = n_side + 1L)
    cuty <- seq(min(data$y), max(data$y), length.out = n_side + 1L)

    category <- interaction(
      findInterval(data$x, cutx, rightmost.closed = TRUE),
      findInterval(data$y, cuty, rightmost.closed = TRUE),
      sep = "-"
    )

    return(factor(category))
  }

  if (cblocks == "row") {

    cuty <- seq(min(data$y), max(data$y), length.out = nblocks + 1L)

    category <- findInterval(data$y, cuty, rightmost.closed = TRUE)
    return(factor(category))
  }

  if (cblocks == "col") {

    cutx <- seq(min(data$x), max(data$x), length.out = nblocks + 1L)

    category <- findInterval(data$x, cutx, rightmost.closed = TRUE)
    return(factor(category))
  }
}
