#' Prepare block-pairwise data
#'
#' @param y numeric response vector.
#' @param coords numeric matrix/data.frame of spatial coordinates.
#' @param X optional design matrix. If NULL, an intercept-only model is used.
#' @param cblocks block construction method.
#' @param nblocks number of spatial blocks.
#' @param fweight pair construction method: "nn" or "dist".
#' @param nnk number of nearest neighbours, used if fweight = "nn".
#' @param w distance cut-off, used if fweight = "dist".
#'
#' @return An object of class `"bp_precomp"` containing block indices,
#'   within-block pairs, pairwise distances, block sizes and response values.
#' @export
prepare_bcl_data <- function(y,
                             coords,
                             X=NULL,
                             cblocks = c("regular", "kmeans", "row", "col"),
                             nblocks = 1L,
                             fweight = c("nn", "dist"),
                             nnk = NULL,
                             w = NULL) {

  cblocks <- match.arg(cblocks)
  fweight <- match.arg(fweight)

  y <- as.numeric(y)
  coords <- as.matrix(coords)

  n <- length(y)

  X <- bp_build_design_matrix(X, n)

  if (length(y) < 2L) {
    stop("y must contain at least two observations.", call. = FALSE)
  }

  if (any(!is.finite(y))) {
    stop("y contains non-finite values.", call. = FALSE)
  }

  if (nrow(coords) != n) {
    stop("nrow(coords) must be equal to length(y).", call. = FALSE)
  }

  sblocks <- set_idblocks_pairslist(
    coords = coords,
    cblocks = cblocks,
    nblocks = nblocks,
    fweight = fweight,
    nnk = nnk,
    w = w
  )

  out <- list(
    id_blocks = sblocks$id_blocks,
    valid_pairs = lapply(sblocks$pairs_list, `[[`, "valid_pairs"),
    dist_pairs = lapply(sblocks$pairs_list, `[[`, "dist_pairs"),
    yk = lapply(sblocks$id_blocks, function(idx) y[idx]),
    coordsk = lapply(sblocks$id_blocks, function(idx) coords[idx, , drop = FALSE]),
    Xk = lapply(sblocks$id_blocks, function(idx) X[idx, , drop = FALSE]),
    groups = sblocks$groups,
    nblocks = length(sblocks$id_blocks),
    block_sizes = lengths(sblocks$id_blocks),
    n_pairs = lengths(lapply(sblocks$pairs_list, `[[`, "dist_pairs")),
    cblocks = cblocks,
    fweight = fweight,
    nnk = nnk,
    w = w
  )

  class(out) <- "bp_precomp"
  out
}


prepare_bcl_ptr <- function(precomp) {
  prepare_bcl_storage(
    y_list_R = precomp$yk,
    X_list_R = precomp$Xk,
    valid_pairs_list_R = precomp$valid_pairs,
    dist_pairs_list_R = precomp$dist_pairs
  )
}


set_idblocks_pairslist <- function(coords,
                                   cblocks = c("regular", "kmeans", "row", "col"),
                                   nblocks = 1L,
                                   fweight = c("nn", "dist"),
                                   nnk = NULL,
                                   w = NULL) {

  cblocks <- match.arg(cblocks)
  fweight <- match.arg(fweight)

  groups <- make_block(coords = coords, nblocks = nblocks, cblocks = cblocks)
  id_blocks <- split(seq_len(nrow(coords)), groups)

  pairs_list <- lapply(id_blocks, function(id_k) {
    coords_k <- as.matrix(coords[id_k, , drop = FALSE])
    make_pairs(coords_k, nnk = nnk, w = w, fweight = fweight)
  })

  list(
    groups = groups,
    id_blocks = id_blocks,
    pairs_list = pairs_list
  )
}
