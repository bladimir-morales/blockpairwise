#' Creating pairs
#'
#' @param coords coordinates of the spatial points with columns x, y.
#' @param fweight pair construction method: "nn" or "dist".
#' @param nnk number of nearest neighbours, used if fweight = "nn".
#' @param w distance cut-off, used if fweight = "dist".
#'
#' @return Object of class list with pairs information.
#' @export
make_pairs <- function(coords,
                       fweight = c("nn", "dist"),
                       nnk = NULL,
                       w = NULL) {

  fweight <- match.arg(fweight)

  coords <- as.matrix(coords)
  if (ncol(coords) != 2L) {
    stop("coords must have exactly two columns.", call. = FALSE)
  }

  n <- nrow(coords)

  if (n < 2L) {
    return(list(
      valid_pairs = matrix(integer(0), ncol = 2L),
      dist_pairs = numeric(0)
    ))
  }

  if (fweight == "dist") {

    if (is.null(w) || !is.finite(w) || w <= 0) {
      stop("w must be positive and finite when fweight = 'dist'.", call. = FALSE)
    }

    D <- as.matrix(stats::dist(coords))
    idx <- which(upper.tri(D) & D < w, arr.ind = TRUE)

    if (nrow(idx) == 0L) {
      return(list(
        valid_pairs = matrix(integer(0), ncol = 2L),
        dist_pairs = numeric(0)
      ))
    }

    list(
      valid_pairs = matrix(idx, ncol = 2L),
      dist_pairs = D[idx]
    )

  } else {

    if (is.null(nnk) || !is.finite(nnk) || nnk < 1L) {
      stop("nnk must be a positive integer when fweight = 'nn'.", call. = FALSE)
    }

    nnk <- as.integer(nnk)
    nnk <- min(nnk, n - 1L)

    nn <- nabor::knn(coords, k = nnk + 1L)
    assemble_pairs(nn$nn.idx, nn$nn.dists, nnk)
  }
}
