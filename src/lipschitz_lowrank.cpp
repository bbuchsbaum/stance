#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Estimate Lipschitz Constant for Low-rank W
//' 
//' Efficiently estimates the Lipschitz constant when W is stored in low-rank form
//' as W = U * diag(S) * V'. Avoids forming the full W matrix.
//' 
//' @param U Left singular vectors (V x r)
//' @param S Singular values (length r)
//' @param V Right singular vectors (K x r)
//' @param hrf_kernel HRF kernel vector
//' @param T Number of time points
//' @param max_iter Maximum iterations for power method
//' @param tol Convergence tolerance
//' 
//' @return Estimated Lipschitz constant
//' 
//' @export
// [[Rcpp::export]]
double estimate_lipschitz_lowrank_rcpp(const arma::mat& U, 
                                      const arma::vec& S,
                                      const arma::mat& V,
                                      const arma::vec& hrf_kernel,
                                      int T,
                                      int max_iter = 30,
                                      double tol = 1e-6) {
  
  // For low-rank W = U * diag(S) * V'
  // W'W = V * diag(S^2) * V'
  // We need the largest eigenvalue of W'W
  
  // Since V * diag(S^2) * V' is already in factored form,
  // the eigenvalues are just S^2 (if V has orthonormal columns)
  // So lambda_max(W'W) = max(S^2)
  
  double lambda_max_WtW = max(square(S));
  
  // Compute ||hrf||_2^2
  double hrf_norm_sq = dot(hrf_kernel, hrf_kernel);
  
  // Return L = ||H||_2^2 * lambda_max(W'W)
  return hrf_norm_sq * lambda_max_WtW;
}

//' Compute W'Y Efficiently Using Low-rank Form
//' 
//' Computes W'Y where W = U * diag(S) * V' without forming full W.
//' 
//' @param U Left singular vectors (V x r)
//' @param S Singular values (length r)
//' @param V Right singular vectors (K x r)
//' @param Y Data matrix (V x T)
//' 
//' @return W'Y matrix (K x T)
//' 
//' @export
// [[Rcpp::export]]
arma::mat compute_WtY_lowrank_rcpp(const arma::mat& U,
                                   const arma::vec& S,
                                   const arma::mat& V,
                                   const arma::mat& Y) {
  // W'Y = V * diag(S) * U' * Y
  // Compute U'Y first (r x T)
  arma::mat UtY = U.t() * Y;
  
  // Scale by S
  for (int i = 0; i < S.n_elem; i++) {
    UtY.row(i) *= S(i);
  }
  
  // Return V * (S * U'Y)
  return V * UtY;
}

//' Compute W'W Efficiently Using Low-rank Form
//' 
//' Computes W'W where W = U * diag(S) * V' without forming full W.
//' 
//' @param V Right singular vectors (K x r)
//' @param S Singular values (length r)
//' 
//' @return W'W matrix (K x K)
//' 
//' @export
// [[Rcpp::export]]
arma::mat compute_WtW_lowrank_rcpp(const arma::mat& V,
                                   const arma::vec& S) {
  // W'W = V * diag(S^2) * V'
  arma::vec S_sq = square(S);
  
  // Scale columns of V by sqrt(S^2) = S
  arma::mat V_scaled = V;
  for (int i = 0; i < S.n_elem; i++) {
    V_scaled.col(i) *= S(i);
  }
  
  // Return V_scaled * V'
  return V_scaled * V.t();
}