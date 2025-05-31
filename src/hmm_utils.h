#ifndef STANCE_HMM_UTILS_H
#define STANCE_HMM_UTILS_H

#include <RcppArmadillo.h>

inline void validate_forward_inputs(const arma::mat& log_lik,
                                    const arma::mat& Pi,
                                    const arma::vec& pi0) {
  if (log_lik.is_empty() || Pi.is_empty() || pi0.is_empty()) {
    Rcpp::stop("Inputs cannot be empty");
  }
  int K = log_lik.n_rows;
  if (Pi.n_rows != K || Pi.n_cols != K) {
    Rcpp::stop("Pi must be K x K");
  }
  if (pi0.n_elem != (unsigned)K) {
    Rcpp::stop("pi0 must have length K");
  }
}

inline void validate_backward_inputs(const arma::mat& log_lik,
                                     const arma::mat& Pi,
                                     const arma::vec& c_scale) {
  if (log_lik.is_empty() || Pi.is_empty() || c_scale.is_empty()) {
    Rcpp::stop("Inputs cannot be empty");
  }
  int K = log_lik.n_rows;
  int T_len = log_lik.n_cols;
  if (Pi.n_rows != K || Pi.n_cols != K) {
    Rcpp::stop("Pi must be K x K");
  }
  if (c_scale.n_elem != (unsigned)T_len) {
    Rcpp::stop("c scale length mismatch");
  }
}

#endif // STANCE_HMM_UTILS_H
