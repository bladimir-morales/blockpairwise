#ifndef PARALLEL_BCL_H
#define PARALLEL_BCL_H

#include <RcppParallel.h>
#include <Rcpp.h>

#include <cstddef>
#include <vector>
#include <thread>
#include <algorithm>
#include <cmath>

#include "bcl_storage.h"


inline int bcl_choose_nthreads(int nthreads) {

  if (nthreads > 0) {
    return std::max(1, nthreads);
  }

  unsigned int detected = std::thread::hardware_concurrency();

  if (detected <= 1U) {
    return 1;
  }

  return std::max(1, static_cast<int>(detected) - 1);
}


inline int bcl_count_pairs(const BCLStorage& storage) {

  int total = 0;

  for (int b = 0; b < storage.nb; ++b) {
    total += static_cast<int>(storage.pair_i[b].size());
  }

  return total;
}


template <class PairContribution>
inline void bcl_validate_storage_contrib(const BCLStorage& storage,
                                         const PairContribution& contrib) {

  if (storage.nb <= 0) {
    Rcpp::stop("storage must contain at least one block.");
  }

  for (int b = 0; b < storage.nb; ++b) {

    const int p = storage.n_cov[b];

    if (p <= 0) {
      Rcpp::stop("ncol(X) must be positive.");
    }

    if (static_cast<int>(contrib.beta.size()) != p) {
      Rcpp::stop("length(beta) must match ncol(X).");
    }

    if (storage.pair_i[b].size() != storage.pair_j[b].size()) {
      Rcpp::stop("pair_i and pair_j must have the same length.");
    }

    if (storage.pair_i[b].size() != storage.dist_pairs_list[b].size()) {
      Rcpp::stop("pair and distance vectors must have the same length.");
    }
  }
}


template <class PairContribution>
inline double eval_block_pairwise_sum(const BCLStorage& storage,
                                      const PairContribution& contrib,
                                      const std::size_t b) {

  const std::vector<double>& yk = storage.y_list[b];
  const std::vector<double>& Xk = storage.X_list[b];

  const int p = storage.n_cov[b];

  const std::vector<int>& pi = storage.pair_i[b];
  const std::vector<int>& pj = storage.pair_j[b];
  const std::vector<double>& dp = storage.dist_pairs_list[b];

  const std::size_t n_pairs = pi.size();

  if (n_pairs == 0) {
    return 0.0;
  }

  const std::size_t nk = yk.size();

  std::vector<double> ek(nk);

  for (std::size_t i = 0; i < nk; ++i) {

    const double* Xi = &Xk[i * static_cast<std::size_t>(p)];

    double mui = 0.0;

    for (int q = 0; q < p; ++q) {
      mui += Xi[q] * contrib.beta[static_cast<std::size_t>(q)];
    }

    ek[i] = yk[i] - mui;
  }

  double local = 0.0;

  for (std::size_t l = 0; l < n_pairs; ++l) {

    const int i1 = pi[l];
    const int i2 = pj[l];

    local += contrib(
      ek[static_cast<std::size_t>(i1)],
        ek[static_cast<std::size_t>(i2)],
          dp[l]
    );
  }

  return local;
}


template <class PairContribution>
struct BlockPairwiseWorker : public RcppParallel::Worker {

  const BCLStorage& storage;
  const PairContribution contrib;
  double total;

  BlockPairwiseWorker(const BCLStorage& storage_,
                      const PairContribution& contrib_)
    : storage(storage_),
      contrib(contrib_),
      total(0.0) {}

  BlockPairwiseWorker(const BlockPairwiseWorker& other,
                      RcppParallel::Split)
    : storage(other.storage),
      contrib(other.contrib),
      total(0.0) {}

  void operator()(std::size_t begin,
                std::size_t end) {

    double local = 0.0;

    for (std::size_t b = begin; b < end; ++b) {
      local += eval_block_pairwise_sum(storage, contrib, b);
    }

    total += local;
  }

  void join(const BlockPairwiseWorker& rhs) {
    total += rhs.total;
  }
};


template <class PairContribution>
double sequential_block_pairwise_sum(const BCLStorage& storage,
                                     const PairContribution& contrib) {

  double total = 0.0;

  for (std::size_t b = 0; b < static_cast<std::size_t>(storage.nb); ++b) {
    total += eval_block_pairwise_sum(storage, contrib, b);
  }

  return total;
}


template <class PairContribution>
double parallel_block_pairwise_sum(const BCLStorage& storage,
                                   const PairContribution& contrib,
                                   int nthreads = 0,
                                   int min_pairs_parallel = 20000) {

  bcl_validate_storage_contrib(storage, contrib);

  const int nt = bcl_choose_nthreads(nthreads);
  const int total_pairs = bcl_count_pairs(storage);

  const bool use_parallel =
    nt > 1 &&
    storage.nb > 1 &&
    total_pairs >= min_pairs_parallel;

  if (!use_parallel) {
    return sequential_block_pairwise_sum(storage, contrib);
  }

  BlockPairwiseWorker<PairContribution> worker(storage, contrib);

  RcppParallel::parallelReduce(
    0,
    static_cast<std::size_t>(storage.nb),
    worker,
    1,
    nt
  );

  return worker.total;
}

#endif
