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


inline int bcl_zi_count_pairs(const BCLZIStorage& storage) {

  int total = 0;

  for (int b = 0; b < storage.nb; ++b) {
    total += static_cast<int>(storage.pair_i[b].size());
  }

  return total;
}


template <class PairZIContribution>
inline void bcl_zi_validate_storage_contrib(const BCLZIStorage& storage,
                                        const PairZIContribution& contrib) {

  if (storage.nb <= 0) {
    Rcpp::stop("ZI storage must contain at least one block.");
  }

  for (int b = 0; b < storage.nb; ++b) {

    if (storage.n_obs[b] <= 0) {
      Rcpp::stop("Each ZI block must contain at least one observation.");
    }

    if (storage.pair_i[b].size() != storage.pair_j[b].size()) {
      Rcpp::stop("pair_i and pair_j must have the same length.");
    }

    if (storage.pair_i[b].size() != storage.dist_pairs_list[b].size()) {
      Rcpp::stop("pair and distance vectors must have the same length.");
    }

    if (contrib.has_X &&
        static_cast<int>(contrib.beta.size()) != storage.p_mu[b]) {
      Rcpp::stop("length(beta) must match ncol(X).");
    }

    if (contrib.has_W &&
        static_cast<int>(contrib.alpha.size()) != storage.p_pi[b]) {
      Rcpp::stop("length(alpha) must match ncol(W).");
    }
  }
}


template <class PairZIContribution>
inline double eval_block_pairwise_zi_sum(const BCLZIStorage& storage,
                                const PairZIContribution& contrib,
                                const std::size_t b) {

  const std::vector<double>& yk = storage.y_list[b];

  const std::vector<double>& Xk = storage.X_list[b];
  const std::vector<double>& Wk = storage.W_list[b];

  const std::vector<int>& pi = storage.pair_i[b];
  const std::vector<int>& pj = storage.pair_j[b];
  const std::vector<double>& dp = storage.dist_pairs_list[b];

  const int n = storage.n_obs[b];
  const int p_mu = storage.p_mu[b];
  const int p_pi = storage.p_pi[b];

  const std::size_t n_pairs = pi.size();

  if (n_pairs == 0) {
    return 0.0;
  }

  std::vector<double> mu(n);
  std::vector<double> prob0(n);

  if (contrib.has_X) {

    for (int i = 0; i < n; ++i) {

      const double* Xi = &Xk[static_cast<std::size_t>(i) *
        static_cast<std::size_t>(p_mu)];

      double eta = 0.0;

      for (int q = 0; q < p_mu; ++q) {
        eta += Xi[q] * contrib.beta[static_cast<std::size_t>(q)];
      }

      mu[i] = std::exp(eta);

      if (!std::isfinite(mu[i]) || mu[i] <= 0.0) {
        return NAN;
      }
    }

  } else {

    const double mu0 = std::exp(contrib.beta[0]);

    if (!std::isfinite(mu0) || mu0 <= 0.0) {
      return NAN;
    }

    for (int i = 0; i < n; ++i) {
      mu[i] = mu0;
    }
  }

  if (contrib.has_W) {

    for (int i = 0; i < n; ++i) {

      const double* Wi = &Wk[static_cast<std::size_t>(i) *
        static_cast<std::size_t>(p_pi)];

      double eta = 0.0;

      for (int q = 0; q < p_pi; ++q) {
        eta += Wi[q] * contrib.alpha[static_cast<std::size_t>(q)];
      }

      prob0[i] = contrib.inv_logit(eta);
    }

  } else {

    const double prob00 = contrib.inv_logit(contrib.alpha[0]);

    for (int i = 0; i < n; ++i) {
      prob0[i] = prob00;
    }
  }

  double local = 0.0;

  for (std::size_t l = 0; l < n_pairs; ++l) {

    const int i = pi[l];
    const int j = pj[l];

    local += contrib(
      yk[static_cast<std::size_t>(i)],
        yk[static_cast<std::size_t>(j)],
          prob0[static_cast<std::size_t>(i)],
               prob0[static_cast<std::size_t>(j)],
                    mu[static_cast<std::size_t>(i)],
                      mu[static_cast<std::size_t>(j)],
                        dp[l]
    );
  }

  return local;
}


template <class PairZIContribution>
struct BlockPairwiseZIWorker : public RcppParallel::Worker {

  const BCLZIStorage& storage;
  const PairZIContribution contrib;
  double total;

  BlockPairwiseZIWorker(const BCLZIStorage& storage_,
                const PairZIContribution& contrib_)
    : storage(storage_),
      contrib(contrib_),
      total(0.0) {}

  BlockPairwiseZIWorker(const BlockPairwiseZIWorker& other,
                RcppParallel::Split)
    : storage(other.storage),
      contrib(other.contrib),
      total(0.0) {}

  void operator()(std::size_t begin,
                std::size_t end) {

    double local = 0.0;

    for (std::size_t b = begin; b < end; ++b) {
      local += eval_block_pairwise_zi_sum(storage, contrib, b);
    }

    total += local;
  }

  void join(const BlockPairwiseZIWorker& rhs) {
    total += rhs.total;
  }
};


template <class PairZIContribution>
double sequential_zi_block_sum(const BCLZIStorage& storage,
                               const PairZIContribution& contrib) {

  double total = 0.0;

  for (std::size_t b = 0; b < static_cast<std::size_t>(storage.nb); ++b) {
    total += eval_block_pairwise_zi_sum(storage, contrib, b);
  }

  return total;
}

template <class PairZIContribution>
double parallel_zi_block_sum(const BCLZIStorage& storage,
                             const PairZIContribution& contrib,
                             int nthreads = 0,
                             int min_pairs_parallel = 20000) {

  bcl_zi_validate_storage_contrib(storage, contrib);

  const int nt = bcl_zi_choose_nthreads(nthreads);
  const int total_pairs = bcl_zi_count_pairs(storage);

  const bool use_parallel =
    nt > 1 &&
    storage.nb > 1 &&
    total_pairs >= min_pairs_parallel;

  if (!use_parallel) {
    return sequential_zi_block_sum(storage, contrib);
  }

  BlockPairwiseZIWorker<PairZIContribution> worker(storage, contrib);

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
