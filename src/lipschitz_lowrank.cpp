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
  // Input validation
  if (U.is_empty() || S.is_empty() || V.is_empty() || Y.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }

  // Dimension checks
  if (U.n_cols != S.n_elem) {
    stop("Dimension mismatch: U.n_cols (%d) != length(S) (%d)",
         U.n_cols, S.n_elem);
  }
  if (V.n_cols != S.n_elem) {
    stop("Dimension mismatch: V.n_cols (%d) != length(S) (%d)",
         V.n_cols, S.n_elem);
  }
  if (Y.n_rows != U.n_rows) {
    stop("Dimension mismatch: Y.n_rows (%d) != U.n_rows (%d)",
         Y.n_rows, U.n_rows);
  }

  // W'Y = V * diag(S) * U' * Y
  // Compute U'Y first (r x T)
  arma::mat UtY = U.t() * Y;

  // Scale rows of U'Y by singular values
  UtY.each_col() %= S;

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
  // Input validation
  if (V.is_empty() || S.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }

  // Dimension check
  if (V.n_cols != S.n_elem) {
    stop("Dimension mismatch: V.n_cols (%d) != length(S) (%d)",
         V.n_cols, S.n_elem);
  }

  // Scale columns of V by singular values
  arma::mat V_scaled = V;
  V_scaled.each_row() %= S.t();

  // Return V_scaled * V_scaled'
  return V_scaled * V_scaled.t();
}
