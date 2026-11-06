#ifndef BCL_ZI_STORAGE_H
#define BCL_ZI_STORAGE_H

#include <Rcpp.h>
#include <vector>
#include <utility>

struct BCLZIStorage {

  std::vector<std::vector<double>> y_list;

  std::vector<int> n_obs;

  std::vector<std::vector<int>> pair_i;
  std::vector<std::vector<int>> pair_j;
  std::vector<std::vector<double>> dist_pairs_list;

  int nb;

  BCLZIStorage(Rcpp::List y_list_R,
               Rcpp::List valid_pairs_list_R,
               Rcpp::List dist_pairs_list_R) {

    nb = y_list_R.size();

    if (valid_pairs_list_R.size() != nb) {
      Rcpp::stop("valid_pairs_list_R must have the same length as y_list_R.");
    }

    if (dist_pairs_list_R.size() != nb) {
      Rcpp::stop("dist_pairs_list_R must have the same length as y_list_R.");
    }

    y_list.reserve(nb);
    n_obs.reserve(nb);
    pair_i.reserve(nb);
    pair_j.reserve(nb);
    dist_pairs_list.reserve(nb);

    for (int b = 0; b < nb; ++b) {

      Rcpp::NumericVector yb =
        Rcpp::as<Rcpp::NumericVector>(y_list_R[b]);

      Rcpp::IntegerMatrix pairs_b =
        Rcpp::as<Rcpp::IntegerMatrix>(valid_pairs_list_R[b]);

      Rcpp::NumericVector dist_b =
        Rcpp::as<Rcpp::NumericVector>(dist_pairs_list_R[b]);

      const int nb_obs = yb.size();
      const int n_pairs = pairs_b.nrow();

      if (pairs_b.ncol() != 2) {
        Rcpp::stop("Each valid_pairs matrix must have exactly 2 columns.");
      }

      if (dist_b.size() != n_pairs) {
        Rcpp::stop("dist_pairs length must match number of valid pairs.");
      }

      std::vector<double> y_cpp(yb.begin(), yb.end());
      std::vector<int> pi_cpp(n_pairs);
      std::vector<int> pj_cpp(n_pairs);
      std::vector<double> dist_cpp(dist_b.begin(), dist_b.end());

      for (int l = 0; l < n_pairs; ++l) {

        int i1 = pairs_b(l, 0) - 1;
        int i2 = pairs_b(l, 1) - 1;

        if (i1 < 0 || i1 >= nb_obs ||
            i2 < 0 || i2 >= nb_obs) {
          Rcpp::stop("Pair index out of bounds.");
        }

        pi_cpp[l] = i1;
        pj_cpp[l] = i2;
      }

      y_list.push_back(std::move(y_cpp));
      n_obs.push_back(nb_obs);
      pair_i.push_back(std::move(pi_cpp));
      pair_j.push_back(std::move(pj_cpp));
      dist_pairs_list.push_back(std::move(dist_cpp));
    }
  }
};


struct BCLZIMarginals {

  std::vector<std::vector<double>> a;
  std::vector<std::vector<double>> b;
  std::vector<std::vector<double>> log1m_pi_h;

  BCLZIMarginals(Rcpp::List marg_blocks,
                 const BCLZIStorage& storage) {

    if (marg_blocks.size() != storage.nb) {
      Rcpp::stop("Number of marginal blocks does not match storage.");
    }

    a.resize(storage.nb);
    b.resize(storage.nb);
    log1m_pi_h.resize(storage.nb);

    for (int k = 0; k < storage.nb; ++k) {

      Rcpp::List marg = marg_blocks[k];

      a[k] = Rcpp::as<std::vector<double>>(marg["a"]);
      b[k] = Rcpp::as<std::vector<double>>(marg["b"]);
      log1m_pi_h[k] =
        Rcpp::as<std::vector<double>>(marg["log1m_pi_h"]);

      const int n = storage.n_obs[k];

      if ((a[k].size() != 1 && a[k].size() != n) ||
          b[k].size() != n ||
          log1m_pi_h[k].size() != n) {
        Rcpp::stop("Invalid marginal dimensions in block.");
      }
    }
  }
};

typedef Rcpp::XPtr<BCLZIStorage> BCLZIStoragePtr;

#endif
