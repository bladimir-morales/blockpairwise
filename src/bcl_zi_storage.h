#ifndef BCL_ZI_STORAGE_H
#define BCL_ZI_STORAGE_H

#include <Rcpp.h>
#include <vector>
#include <utility>

struct BCLZIStorage {

  std::vector<std::vector<double>> y_list;

  std::vector<std::vector<double>> X_list;
  std::vector<std::vector<double>> W_list;

  std::vector<int> n_obs;
  std::vector<int> p_mu;
  std::vector<int> p_pi;

  std::vector<std::vector<int>> pair_i;
  std::vector<std::vector<int>> pair_j;
  std::vector<std::vector<double>> dist_pairs_list;

  int nb;

  BCLZIStorage(Rcpp::List y_list_R,
               Rcpp::List valid_pairs_list_R,
               Rcpp::List dist_pairs_list_R,
               Rcpp::Nullable<Rcpp::List> X_list_R = R_NilValue,
               Rcpp::Nullable<Rcpp::List> W_list_R = R_NilValue,
               bool has_X = false,
               bool has_W = false) {

    nb = y_list_R.size();

    if (valid_pairs_list_R.size() != nb) {
      Rcpp::stop("valid_pairs_list_R must have the same length as y_list_R.");
    }

    if (dist_pairs_list_R.size() != nb) {
      Rcpp::stop("dist_pairs_list_R must have the same length as y_list_R.");
    }

    Rcpp::List Xs;
    Rcpp::List Ws;

    if (has_X) {
      Xs = Rcpp::List(X_list_R);

      if (Xs.size() != nb) {
        Rcpp::stop("X_list_R must have the same length as y_list_R.");
      }
    }

    if (has_W) {
      Ws = Rcpp::List(W_list_R);

      if (Ws.size() != nb) {
        Rcpp::stop("W_list_R must have the same length as y_list_R.");
      }
    }

    y_list.reserve(nb);
    X_list.reserve(nb);
    W_list.reserve(nb);

    n_obs.reserve(nb);
    p_mu.reserve(nb);
    p_pi.reserve(nb);

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

      int p_mu_b = 1;
      int p_pi_b = 1;

      std::vector<double> X_cpp;
      std::vector<double> W_cpp;

      if (has_X) {

        Rcpp::NumericMatrix Xb =
          Rcpp::as<Rcpp::NumericMatrix>(Xs[b]);

        if (Xb.nrow() != nb_obs) {
          Rcpp::stop("Each X matrix must have the same number of rows as y.");
        }

        p_mu_b = Xb.ncol();

        if (p_mu_b <= 0) {
          Rcpp::stop("X matrices must have at least one column.");
        }

        X_cpp.resize(nb_obs * p_mu_b);

        for (int i = 0; i < nb_obs; ++i) {
          for (int q = 0; q < p_mu_b; ++q) {
            X_cpp[i * p_mu_b + q] = Xb(i, q);
          }
        }
      }

      if (has_W) {

        Rcpp::NumericMatrix Wb =
          Rcpp::as<Rcpp::NumericMatrix>(Ws[b]);

        if (Wb.nrow() != nb_obs) {
          Rcpp::stop("Each W matrix must have the same number of rows as y.");
        }

        p_pi_b = Wb.ncol();

        if (p_pi_b <= 0) {
          Rcpp::stop("W matrices must have at least one column.");
        }

        W_cpp.resize(nb_obs * p_pi_b);

        for (int i = 0; i < nb_obs; ++i) {
          for (int q = 0; q < p_pi_b; ++q) {
            W_cpp[i * p_pi_b + q] = Wb(i, q);
          }
        }
      }

      std::vector<int> pi_cpp(n_pairs);
      std::vector<int> pj_cpp(n_pairs);
      std::vector<double> dist_cpp(dist_b.begin(), dist_b.end());

      for (int l = 0; l < n_pairs; ++l) {

        int i1 = pairs_b(l, 0) - 1;
        int i2 = pairs_b(l, 1) - 1;

        if (i1 < 0 || i1 >= nb_obs || i2 < 0 || i2 >= nb_obs) {
          Rcpp::stop("Pair index out of bounds.");
        }

        pi_cpp[l] = i1;
        pj_cpp[l] = i2;
      }

      y_list.push_back(std::move(y_cpp));
      X_list.push_back(std::move(X_cpp));
      W_list.push_back(std::move(W_cpp));

      n_obs.push_back(nb_obs);
      p_mu.push_back(p_mu_b);
      p_pi.push_back(p_pi_b);

      pair_i.push_back(std::move(pi_cpp));
      pair_j.push_back(std::move(pj_cpp));
      dist_pairs_list.push_back(std::move(dist_cpp));
    }
  }
};

typedef Rcpp::XPtr<BCLZIStorage> BCLZIStoragePtr;

#endif
