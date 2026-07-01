#ifndef PARALLEL_BCL_H
#define PARALLEL_BCL_H

#include <RcppParallel.h>
#include <Rcpp.h>
#include <cstddef>
#include <vector>

#include "bcl_storage.h"


// ============================================================
// 1. Worker general con covariables
// ============================================================

template <class PairContribution>
struct BlockPairwiseWorker : public RcppParallel::Worker {

  const BCLStorage& storage;
  const PairContribution contrib;
  double total;

  BlockPairwiseWorker(const BCLStorage& storage_,
                      const PairContribution& contrib_)
    : storage(storage_), contrib(contrib_), total(0.0) {}

  BlockPairwiseWorker(const BlockPairwiseWorker& other,
                      RcppParallel::Split)
    : storage(other.storage), contrib(other.contrib), total(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {

    double local = 0.0;

    for (std::size_t b = begin; b < end; ++b) {

      const std::vector<double>& yk = storage.y_list[b];
      const std::vector<double>& Xk = storage.X_list[b];
      const int p = storage.n_cov[b];

      const std::vector<int>& pi = storage.pair_i[b];
      const std::vector<int>& pj = storage.pair_j[b];
      const std::vector<double>& dp = storage.dist_pairs_list[b];

      const std::size_t n_pairs = pi.size();
      if (n_pairs == 0) continue;

      const std::size_t nk = yk.size();

      if (p <= 0) {
        Rcpp::stop("ncol(X) must be positive.");
      }

      if (static_cast<int>(contrib.beta.size()) != p) {
        Rcpp::stop("length(beta) must match ncol(X).");
      }

      std::vector<double> ek(nk);

      for (std::size_t i = 0; i < nk; ++i) {

        const double* Xi = &Xk[i * static_cast<std::size_t>(p)];

        double mui = 0.0;

        for (int q = 0; q < p; ++q) {
          mui += Xi[q] * contrib.beta[static_cast<std::size_t>(q)];
        }

        ek[i] = yk[i] - mui;
      }

      for (std::size_t i = 0; i < n_pairs; ++i) {

        const int i1 = pi[i];
        const int i2 = pj[i];

        local += contrib(
          ek[static_cast<std::size_t>(i1)],
            ek[static_cast<std::size_t>(i2)],
              dp[i]
        );
      }
    }

    total += local;
  }

  void join(const BlockPairwiseWorker& rhs) {
    total += rhs.total;
  }
};



template <class PairContribution>
double parallel_block_pairwise_sum(const BCLStorage& storage,
                                   const PairContribution& contrib) {

  BlockPairwiseWorker<PairContribution> worker(storage, contrib);

  RcppParallel::parallelReduce(
    0,
    static_cast<std::size_t>(storage.nb),
    worker
  );

  return worker.total;
}

#endif
