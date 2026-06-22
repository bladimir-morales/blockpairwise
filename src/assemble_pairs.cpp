#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List assemble_pairs(IntegerMatrix nn_idx, NumericMatrix nn_dists, int nnk) {

  if (nnk <= 0) {
    stop("nnk should be positive.");
  }

  if (nn_idx.ncol() < nnk + 1) {
    stop("nn_idx should have at least nnk + 1 columns.");
  }

  if (nn_dists.ncol() < nnk + 1) {
    stop("nn_dists should have at least nnk + 1 columns.");
  }

  int n = nn_idx.nrow();
  int total = n * nnk;

  IntegerMatrix valid_pairs(total, 2);
  NumericVector dist_pairs(total);

  int idx = 0;
  for (int i = 0; i < n; i++) {
    int self = nn_idx(i, 0);
    for (int j = 1; j <= nnk; j++) {
      valid_pairs(idx, 0) = self;
      valid_pairs(idx, 1) = nn_idx(i, j);
      dist_pairs[idx] = nn_dists(i, j);
      idx++;
    }
  }
  return List::create(
    _["valid_pairs"] = valid_pairs,
    _["dist_pairs"]  = dist_pairs
  );
}

