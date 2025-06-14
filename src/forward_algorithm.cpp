#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

// Forward algorithm used by the ContinuousBayesianDecoder.  The
// inputs are typically generated by the Rcpp log-likelihood kernel and
// work for matrices obtained from neuroim2 objects.
#include "hmm_utils.h"

using namespace Rcpp;
using namespace arma;

//' Forward Algorithm for HMM (CBD)
//'
//' Computes the forward probabilities (alpha) for a sequence of
//' log-likelihoods under an HMM with transition matrix `Pi` and
//' initial distribution `pi0`. Scaling factors are used for numerical
//' stability and the overall log-likelihood of the sequence is
//' returned.  It mirrors the R implementation but benefits from
//' Armadillo matrix operations.
//'
//' @param log_lik K x T matrix of log-likelihoods
//' @param Pi      K x K state transition matrix
//' @param pi0     Vector of length K of initial state probabilities
//'
//' @return List with `alpha` matrix (K x T), `log_likelihood`, and
//'   scaling vector `c` used during the recursion.
//'
//' @keywords internal
// [[Rcpp::export]]
List forward_pass_rcpp(const arma::mat& log_lik,
                       const arma::mat& Pi,
                       const arma::vec& pi0) {
  validate_forward_inputs(log_lik, Pi, pi0);

  int K = log_lik.n_rows;
  int T_len = log_lik.n_cols;

  arma::mat alpha(K, T_len, fill::zeros);
  arma::vec c_scale(T_len, fill::ones);

  // initialize
  arma::vec alpha_col = pi0 % exp(log_lik.col(0));
  double denom = accu(alpha_col);
  if (denom <= 0.0) denom = 1e-12;
  c_scale(0) = 1.0 / denom;
  alpha_col *= c_scale(0);
  alpha.col(0) = alpha_col;

  for (int t = 1; t < T_len; ++t) {
    alpha_col = (Pi.t() * alpha.col(t - 1));
    alpha_col %= exp(log_lik.col(t));
    denom = accu(alpha_col);
    if (denom <= 0.0) denom = 1e-12;
    c_scale(t) = 1.0 / denom;
    alpha_col *= c_scale(t);
    alpha.col(t) = alpha_col;
  }

  double log_likelihood = -sum(log(c_scale));

  return List::create(Named("alpha") = alpha,
                      Named("log_likelihood") = log_likelihood,
                      Named("c") = c_scale);
}

