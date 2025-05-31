#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Companion to `forward_pass_rcpp`; processes the sequence in
// reverse using Armadillo operations.  Designed for use with data
// converted from neuroim2 objects.
#include "hmm_utils.h"

using namespace Rcpp;
using namespace arma;

//' Backward Algorithm for HMM (CBD)
//'
//' Computes the backward probabilities (beta) for a sequence of
//' log-likelihoods given scaling factors from the forward pass.  The
//' implementation mirrors the pure R version but is significantly
//' faster for large matrices.
//'
//' @param log_lik K x T matrix of log-likelihoods
//' @param Pi      K x K state transition matrix
//' @param c_scale Scaling factors computed during forward pass
//'
//' @return Matrix beta (K x T)
//' @keywords internal
// [[Rcpp::export]]
arma::mat backward_pass_rcpp(const arma::mat& log_lik,
                             const arma::mat& Pi,
                             const arma::vec& c_scale) {
  validate_backward_inputs(log_lik, Pi, c_scale);

  int K = log_lik.n_rows;
  int T_len = log_lik.n_cols;

  arma::mat beta(K, T_len, fill::ones);
  beta.col(T_len - 1).ones();

  if (T_len > 1) {
    for (int t = T_len - 2; t >= 0; --t) {
      arma::vec tmp = beta.col(t + 1) % exp(log_lik.col(t + 1));
      beta.col(t) = Pi * tmp;
      beta.col(t) *= c_scale(t + 1);
    }
  }

  return beta;
}
