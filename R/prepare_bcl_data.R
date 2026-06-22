#' Prepare block-pairwise data
#'
#' @param data data.frame or matrix with columns x, y, z.
#' @param cblocks block construction method.
#' @param nblocks number of spatial blocks.
#' @param fweight pair construction method: "nn" or "dist".
#' @param nnk number of nearest neighbours, used if fweight = "nn".
#' @param w distance cut-off, used if fweight = "dist".
#'
#' @return An object of class `"bp_precomp"` containing block indices,
#'   within-block pairs, pairwise distances, block sizes and response values.
#' @export
prepare_bcl_data <- function(data,
                             cblocks = c("regular", "kmeans", "row", "col"),
                             nblocks = 1L,
                             fweight = c("nn", "dist"),
                             nnk = NULL,
                             w = NULL) {

  data <- validate_spatial_data(data)
  cblocks <- match.arg(cblocks)
  fweight <- match.arg(fweight)

  sblocks <- set_idblocks_pairslist(
    data = data,
    cblocks = cblocks,
    nblocks = nblocks,
    fweight = fweight,
    nnk = nnk,
    w = w
  )

  out <- list(
    data = data,
    id_blocks = sblocks$id_blocks,
    valid_pairs = lapply(sblocks$pairs_list, `[[`, "valid_pairs"),
    dist_pairs = lapply(sblocks$pairs_list, `[[`, "dist_pairs"),
    Zk = lapply(sblocks$id_blocks, function(idx) data$z[idx]),
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
    Z_list_R = precomp$Zk,
    valid_pairs_list_R = precomp$valid_pairs,
    dist_pairs_list_R = precomp$dist_pairs
  )
}


set_idblocks_pairslist <- function(data,
                                   cblocks = c("regular", "kmeans", "row", "col"),
                                   nblocks = 1L,
                                   fweight = c("nn", "dist"),
                                   nnk = NULL,
                                   w = NULL) {

  data <- validate_spatial_data(data)
  cblocks <- match.arg(cblocks)
  fweight <- match.arg(fweight)

  groups <- make_block(data, nblocks = nblocks, cblocks = cblocks)
  id_blocks <- split(seq_len(nrow(data)), groups)

  pairs_list <- lapply(id_blocks, function(id_k) {
    coords_k <- as.matrix(data[id_k, c("x", "y"), drop = FALSE])
    make_pairs(coords_k, nnk = nnk, w = w, fweight = fweight)
  })

  list(
    groups = groups,
    id_blocks = id_blocks,
    pairs_list = pairs_list
  )
}
