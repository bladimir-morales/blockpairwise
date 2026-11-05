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

  const auto& pik = marg.pi[k];
  const auto& ak = marg.a[k];
  const auto& bk = marg.b[k];
  const auto& Lk = marg.log1m_pi_h[k];

  double local = 0.0;

  for (std::size_t l = 0; l < pair_i.size(); ++l) {

    const std::size_t i = static_cast<std::size_t>(pair_i[l]);
    const std::size_t j = static_cast<std::size_t>(pair_j[l]);

    const double pi_i = pik.size() == 1 ? pik[0] : pik[i];
    const double pi_j = pik.size() == 1 ? pik[0] : pik[j];

    const double a_i = ak.size() == 1 ? ak[0] : ak[i];
    const double a_j = ak.size() == 1 ? ak[0] : ak[j];

    local += contrib(
      yk[i], yk[j],
               pi_i, pi_j,
               a_i, a_j,
               bk[i], bk[j],
                        Lk[i], Lk[j],
                                 dist[l]
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
                             int nthreads = 0) {
  bcl_zi_validate_storage_contrib(storage, contrib);

  const int nt_requested = bcl_zi_choose_nthreads(nthreads);
  const int nt = std::min(nt_requested, storage.nb);

  const bool use_parallel =
    nt > 1 &&
    storage.nb > 1;

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
