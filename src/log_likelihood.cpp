#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;

//' Compute Log Likelihoods (CBD)
//'
//' Calculates log-likelihood values for each latent state at each
//' time point. Inputs are already projected into the low-rank
//' space to avoid expensive V x T operations.
//'
//' @param Y_proj r x T matrix of projected data (U'Y)
//' @param Vmat   K x r matrix of state loadings
//' @param hrf_kernel HRF kernel vector
//' @param sigma2 Shared noise variance
//'
//' @return K x T matrix of log-likelihoods
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_log_likelihoods_rcpp(const arma::mat& Y_proj,
                                       const arma::mat& Vmat,
                                       const arma::vec& hrf_kernel,
                                       double sigma2) {
  if (Y_proj.is_empty() || Vmat.is_empty() || hrf_kernel.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }

  int r = Y_proj.n_rows;
  int T_len = Y_proj.n_cols;
  int K = Vmat.n_rows;

  if (Vmat.n_cols != r) {
    stop("Dimension mismatch: Vmat.n_cols (%d) != Y_proj.n_rows (%d)",
         Vmat.n_cols, r);
  }

  sigma2 = std::max(sigma2, 1e-8);
  double const_term = -0.5 * r * std::log(2.0 * M_PI * sigma2);

  // First element of HRF (equivalent to diag convolution in R code)
  double h0 = hrf_kernel.n_elem > 0 ? hrf_kernel(0) : 0.0;
  arma::rowvec h_at_t(T_len, fill::value(h0));

  arma::mat log_lik(K, T_len, fill::zeros);

  for (int k = 0; k < K; ++k) {
    arma::mat mu = Vmat.row(k).t() * h_at_t;
    arma::mat diff = Y_proj - mu;
    arma::rowvec colsum = sum(square(diff), 0);
    log_lik.row(k) = const_term - 0.5 * colsum / sigma2;
  }

  return log_lik;
}

