#ifndef PARALLEL_BCL_ZI_H
#define PARALLEL_BCL_ZI_H

#include <RcppParallel.h>
#include <Rcpp.h>

#include <vector>
#include <thread>
#include <algorithm>
#include <cmath>
#include <cstddef>

#include "bcl_zi_storage.h"


inline int bcl_zi_choose_nthreads(int nthreads) {

  if (nthreads > 0) {
    return std::max(1, nthreads);
  }

  unsigned int detected = std::thread::hardware_concurrency();

  if (detected <= 1U) {
    return 1;
  }

  return std::max(1, static_cast<int>(detected) - 1);
}



// ============================================================
// Block pairwise evaluation with precomputed ZI marginals
// ============================================================

template <class PairZIContribution>
inline double eval_block_pairwise_zi_marginals_sum(
    const BCLZIStorage& storage,
    const BCLZIMarginals& marg,
    const PairZIContribution& contrib,
    const std::size_t k) {

  const auto& yk = storage.y_list[k];
  const auto& pair_i = storage.pair_i[k];
  const auto& pair_j = storage.pair_j[k];
  const auto& dist = storage.dist_pairs_list[k];

  const auto& ak = marg.a[k];
  const auto& bk = marg.b[k];
  const auto& Lk = marg.log1m_pi_h[k];

  double local = 0.0;

  for (std::size_t l = 0; l < pair_i.size(); ++l) {

    const std::size_t i = static_cast<std::size_t>(pair_i[l]);
    const std::size_t j = static_cast<std::size_t>(pair_j[l]);

    const double a_i = ak.size() == 1 ? ak[0] : ak[i];
    const double a_j = ak.size() == 1 ? ak[0] : ak[j];

    local += contrib(yk[i], yk[j], a_i, a_j,
                     bk[i], bk[j], Lk[i], Lk[j], dist[l]
    );
  }

  return local;
}



// ============================================================
// Parallel worker with precomputed ZI marginals
// ============================================================

template <class PairZIContribution>
struct BlockPairwiseZIMarginalsWorker : public RcppParallel::Worker {

  const BCLZIStorage& storage;
  const BCLZIMarginals& marg;
  const PairZIContribution contrib;
  double total;

  BlockPairwiseZIMarginalsWorker(
    const BCLZIStorage& storage_,
    const BCLZIMarginals& marg_,
    const PairZIContribution& contrib_)
    : storage(storage_), marg(marg_), contrib(contrib_), total(0.0) {}

  BlockPairwiseZIMarginalsWorker(
    const BlockPairwiseZIMarginalsWorker& other,
    RcppParallel::Split)
    : storage(other.storage), marg(other.marg),
      contrib(other.contrib), total(0.0) {}

  void operator()(std::size_t begin, std::size_t end) {

    double local = 0.0;

    for (std::size_t k = begin; k < end; ++k) {
      local += eval_block_pairwise_zi_marginals_sum(
        storage, marg, contrib, k
      );
    }

    total += local;
  }

  void join(const BlockPairwiseZIMarginalsWorker& rhs) {
    total += rhs.total;
  }
};


// ============================================================
// Parallel evaluation with precomputed marginals
// ============================================================

template <class PairZIContribution>
double parallel_zi_marginals_block_sum(
    const BCLZIStorage& storage,
    const BCLZIMarginals& marg,
    const PairZIContribution& contrib,
    int nthreads = 0) {

  if (storage.nb <= 0) {
    Rcpp::stop("ZI storage must contain at least one block.");
  }

  const int nt = std::min(
    bcl_zi_choose_nthreads(nthreads), storage.nb
  );

  if (nt <= 1) {
    double total = 0.0;

    for (std::size_t k = 0; k < static_cast<std::size_t>(storage.nb); ++k) {
      total += eval_block_pairwise_zi_marginals_sum(
        storage, marg, contrib, k
      );
    }

    return total;
  }

  BlockPairwiseZIMarginalsWorker<PairZIContribution> worker(
      storage, marg, contrib
  );

  RcppParallel::parallelReduce(
    0, static_cast<std::size_t>(storage.nb),
    worker, 1, nt
  );

  return worker.total;
}
#endif
