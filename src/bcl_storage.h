#ifndef BCL_STORAGE_H
#define BCL_STORAGE_H

#include <Rcpp.h>
#include <vector>

struct BCLStorage {

  std::vector<std::vector<double>> Z_list;

  // X por bloque, guardado en formato fila-major:
  // X_list[b][i * p + q] = X_b(i, q)
  std::vector<std::vector<double>> X_list;

  std::vector<int> n_obs;
  std::vector<int> n_cov;

  std::vector<std::vector<int>> pair_i;
  std::vector<std::vector<int>> pair_j;
  std::vector<std::vector<double>> dist_pairs_list;

  int nb;

  BCLStorage(Rcpp::List Z_list_R,
             Rcpp::List X_list_R,
             Rcpp::List valid_pairs_list_R,
             Rcpp::List dist_pairs_list_R) {

    nb = Z_list_R.size();

    if (X_list_R.size() != nb) {
      Rcpp::stop("X_list_R must have the same length as Z_list_R.");
    }

    if (valid_pairs_list_R.size() != nb) {
      Rcpp::stop("valid_pairs_list_R must have the same length as Z_list_R.");
    }

    if (dist_pairs_list_R.size() != nb) {
      Rcpp::stop("dist_pairs_list_R must have the same length as Z_list_R.");
    }

    Z_list.reserve(nb);
    X_list.reserve(nb);
    n_obs.reserve(nb);
    n_cov.reserve(nb);

    pair_i.reserve(nb);
    pair_j.reserve(nb);
    dist_pairs_list.reserve(nb);

    for (int b = 0; b < nb; ++b) {

      Rcpp::NumericVector yb =
        Rcpp::as<Rcpp::NumericVector>(Z_list_R[b]);

      Rcpp::NumericMatrix Xb =
        Rcpp::as<Rcpp::NumericMatrix>(X_list_R[b]);

      Rcpp::IntegerMatrix pairs_b =
        Rcpp::as<Rcpp::IntegerMatrix>(valid_pairs_list_R[b]);

      Rcpp::NumericVector dist_b =
        Rcpp::as<Rcpp::NumericVector>(dist_pairs_list_R[b]);

      const int nb_obs = yb.size();
      const int p = Xb.ncol();
      const int n_pairs = pairs_b.nrow();

      if (Xb.nrow() != nb_obs) {
        Rcpp::stop("Each X matrix must have the same number of rows as its y vector.");
      }

      if (pairs_b.ncol() != 2) {
        Rcpp::stop("Each valid_pairs matrix must have exactly 2 columns.");
      }

      if (dist_b.size() != n_pairs) {
        Rcpp::stop("dist_pairs length must match number of valid pairs.");
      }

      std::vector<double> y_cpp(yb.begin(), yb.end());

      std::vector<double> X_cpp(nb_obs * p);

      for (int i = 0; i < nb_obs; ++i) {
        for (int q = 0; q < p; ++q) {
          X_cpp[i * p + q] = Xb(i, q);
        }
      }

      std::vector<int> pi_cpp(n_pairs);
      std::vector<int> pj_cpp(n_pairs);
      std::vector<double> dist_cpp(dist_b.begin(), dist_b.end());

      for (int i = 0; i < n_pairs; ++i) {

        int i1 = pairs_b(i, 0) - 1;
        int i2 = pairs_b(i, 1) - 1;

        if (i1 < 0 || i1 >= nb_obs || i2 < 0 || i2 >= nb_obs) {
          Rcpp::stop("Pair index out of bounds.");
        }

        pi_cpp[i] = i1;
        pj_cpp[i] = i2;
      }

      Z_list.push_back(std::move(y_cpp));
      X_list.push_back(std::move(X_cpp));

      n_obs.push_back(nb_obs);
      n_cov.push_back(p);

      pair_i.push_back(std::move(pi_cpp));
      pair_j.push_back(std::move(pj_cpp));
      dist_pairs_list.push_back(std::move(dist_cpp));
    }
  }
};

typedef Rcpp::XPtr<BCLStorage> BCLPtr;

#endif
