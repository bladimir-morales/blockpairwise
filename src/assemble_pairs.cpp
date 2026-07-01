#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List assemble_pairs(IntegerMatrix nn_idx,
                    NumericMatrix nn_dists,
                    int nnk) {

  if (nnk <= 0) {
    stop("nnk should be positive.");
  }

  const int n = nn_idx.nrow();

  if (n <= 0) {
    stop("nn_idx should have at least one row.");
  }

  if (nn_dists.nrow() != n) {
    stop("nn_dists must have the same number of rows as nn_idx.");
  }

  if (nn_idx.ncol() < nnk + 1) {
    stop("nn_idx should have at least nnk + 1 columns.");
  }

  if (nn_dists.ncol() < nnk + 1) {
    stop("nn_dists should have at least nnk + 1 columns.");
  }

  std::vector<int> pair_i;
  std::vector<int> pair_j;
  std::vector<double> pair_d;

  pair_i.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(nnk));
  pair_j.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(nnk));
  pair_d.reserve(static_cast<std::size_t>(n) * static_cast<std::size_t>(nnk));

  for (int i = 0; i < n; ++i) {

    const int self = i + 1;  // R indexing
    int added = 0;

    for (int j = 0; j < nn_idx.ncol(); ++j) {

      if (added >= nnk) {
        break;
      }

      const int neigh = nn_idx(i, j);
      const double dij = nn_dists(i, j);

      if (neigh == NA_INTEGER) {
        continue;
      }

      if (neigh < 1 || neigh > n) {
        stop("nn_idx contains an index outside 1:n.");
      }

      if (neigh == self) {
        continue;
      }

      if (!R_finite(dij)) {
        stop("nn_dists contains non-finite distances.");
      }

      if (dij < 0.0) {
        stop("nn_dists contains negative distances.");
      }

      pair_i.push_back(self);
      pair_j.push_back(neigh);
      pair_d.push_back(dij);

      ++added;
    }
  }

  const int total = static_cast<int>(pair_i.size());

  IntegerMatrix valid_pairs(total, 2);
  NumericVector dist_pairs(total);

  for (int k = 0; k < total; ++k) {
    valid_pairs(k, 0) = pair_i[static_cast<std::size_t>(k)];
    valid_pairs(k, 1) = pair_j[static_cast<std::size_t>(k)];
    dist_pairs[k] = pair_d[static_cast<std::size_t>(k)];
  }

  return List::create(
    _["valid_pairs"] = valid_pairs,
    _["dist_pairs"]  = dist_pairs
  );
}
