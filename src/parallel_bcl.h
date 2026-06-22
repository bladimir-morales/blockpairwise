#ifndef PARALLEL_BCL_H
#define PARALLEL_BCL_H

#include <RcppParallel.h>
#include <cstddef>
#include "bcl_storage.h"

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

      const std::vector<double>& Zk = storage.Z_list[b];
      const std::vector<int>& pi = storage.pair_i[b];
      const std::vector<int>& pj = storage.pair_j[b];
      const std::vector<double>& dp = storage.dist_pairs_list[b];

      const std::size_t n_pairs = pi.size();

      if (n_pairs == 0) continue;

      for (std::size_t i = 0; i < n_pairs; ++i) {

        const int i1 = pi[i];
        const int i2 = pj[i];

        local += contrib(Zk[i1], Zk[i2], dp[i]);
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
