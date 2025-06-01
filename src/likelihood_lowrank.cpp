#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

// Forward declaration of existing row-wise convolution
arma::mat convolve_rows_rcpp(const arma::mat& X, const arma::vec& hrf,
                             int n_threads = 0);

//' Complete low-rank log-likelihood calculation
//'
//' Computes log-likelihood values using the low-rank projection
//' of the data. This function extends the simpler
//' `compute_log_likelihoods_rcpp()` by convolving the state
//' probabilities with the HRF kernel and working entirely in the
//' projected r-dimensional space.
//'
//' @param Y_proj r x T matrix U'Y
//' @param U V x r orthonormal spatial basis
//' @param Vmat K x r matrix of loadings
//' @param H_v V x L_basis HRF coefficients
//' @param hrf_basis L_h x L_basis HRF basis matrix
//' @param S_gamma K x T matrix of state posteriors
//' @param sigma2 noise variance
//'
//' @return K x T matrix of log-likelihoods
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_log_likelihood_lowrank_complete(const arma::mat& Y_proj,
                                                  const arma::mat& U,
                                                  const arma::mat& Vmat,
                                                  const arma::mat& H_v,
                                                  const arma::mat& hrf_basis,
                                                  const arma::mat& S_gamma,
                                                  double sigma2) {
  if (Y_proj.is_empty() || Vmat.is_empty() || H_v.is_empty() ||
      hrf_basis.is_empty() || S_gamma.is_empty()) {
    stop("Inputs cannot be empty");
  }

  int r = Y_proj.n_rows;
  int T_len = Y_proj.n_cols;
  int K = Vmat.n_rows;

  // Average HRF across voxels
  arma::vec hrf_global = hrf_basis * arma::mean(H_v, 0).t();

  // Convolve state probabilities with HRF
  arma::mat S_conv = convolve_rows_rcpp(S_gamma, hrf_global, 0);

  sigma2 = std::max(sigma2, 1e-8);
  double const_term = -0.5 * r * std::log(2.0 * M_PI * sigma2);

  arma::mat log_lik(K, T_len, fill::zeros);

  for (int t = 0; t < T_len; ++t) {
    arma::vec y_t = Y_proj.col(t);
    for (int k = 0; k < K; ++k) {
      arma::vec mu = Vmat.row(k).t() * S_conv(k, t);
      arma::vec diff = y_t - mu;
      double sq_err = arma::dot(diff, diff);
      log_lik(k, t) = const_term - 0.5 * sq_err / sigma2;
    }
  }

  return log_lik;
}

