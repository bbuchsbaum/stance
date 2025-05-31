#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Proximal Operator for Total Variation Penalty (Condat's Algorithm)
//' 
//' Computes the proximal operator for the 1D Total Variation penalty
//' using Condat's direct algorithm. This solves:
//' argmin_z { 0.5 ||z - x||^2 + lambda * TV(z) }
//' where TV(z) = sum_i |z[i] - z[i-1]|
//' 
//' @param x Input vector
//' @param lambda TV regularization parameter (scaled by step size)
//' 
//' @return Solution vector
//' 
//' @references
//' Condat, L. (2013). A direct algorithm for 1-D total variation denoising.
//' IEEE Signal Processing Letters, 20(11), 1054-1057.
//' 
//' @export
// [[Rcpp::export]]
arma::vec prox_tv_condat_1d(const arma::vec& x, double lambda) {
  // Input validation
  if (x.is_empty()) {
    stop("Input vector cannot be empty");
  }
  
  int n = x.n_elem;
  
  if (n <= 1 || lambda <= 0) {
    return x;
  }
  
  // Check for numerical issues
  if (x.has_nan() || x.has_inf()) {
    stop("Input vector contains NaN or Inf values");
  }
  if (!std::isfinite(lambda)) {
    stop("Lambda must be finite");
  }
  
  // Condat's algorithm implementation
  arma::vec z(n);
  
  // Working variables
  double mn = x(0) - lambda;
  double mx = x(0) + lambda;
  int j = 0;
  
  for (int i = 1; i < n; i++) {
    if (x(i) + lambda < mn) {
      // Extend current segment down
      for (int k = j; k <= i - 1; k++) {
        z(k) = mn;
      }
      j = i;
      mn = x(i) - lambda;
      mx = x(i) + lambda;
    } else if (x(i) - lambda > mx) {
      // Extend current segment up
      for (int k = j; k <= i - 1; k++) {
        z(k) = mx;
      }
      j = i;
      mn = x(i) - lambda;
      mx = x(i) + lambda;
    } else {
      // Update bounds
      if (x(i) - lambda > mn) {
        mn = x(i) - lambda;
      }
      if (x(i) + lambda < mx) {
        mx = x(i) + lambda;
      }
    }
  }
  
  // Handle last segment
  double segment_mean = 0;
  int segment_length = n - j;
  if (segment_length <= 0) {
    stop("Invalid segment length in Condat algorithm");
  }
  
  for (int k = j; k < n; k++) {
    segment_mean += x(k);
  }
  segment_mean /= segment_length;
  
  // Clip to feasible range
  if (segment_mean < mn) segment_mean = mn;
  if (segment_mean > mx) segment_mean = mx;
  
  for (int k = j; k < n; k++) {
    z(k) = segment_mean;
  }
  
  return z;
}

//' Proximal Operator for TV Penalty on Matrix Rows
//' 
//' Applies the 1D TV proximal operator to each row of a matrix.
//' This is used in FISTA to enforce temporal smoothness on state activations.
//' Supports OpenMP parallelization when available.
//' 
//' @param X Input matrix (K x T)
//' @param lambda_tv TV regularization parameter
//' @param n_threads Number of threads to use (0 = auto)
//' 
//' @return Matrix with TV-denoised rows
//' 
//' @export
// [[Rcpp::export]]
arma::mat prox_tv_condat_rcpp(const arma::mat& X, double lambda_tv,
                              int n_threads = 0) {
  // Input validation
  if (X.is_empty()) {
    stop("Input matrix cannot be empty");
  }
  if (!std::isfinite(lambda_tv) || lambda_tv < 0) {
    stop("lambda_tv must be non-negative and finite");
  }
  
  int K = X.n_rows;
  int T = X.n_cols;

#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n_threads <= 0 || n_threads > max_threads) {
    n_threads = max_threads;
  }
  omp_set_num_threads(n_threads);
#endif
  
  // Check for numerical issues
  if (X.has_nan() || X.has_inf()) {
    stop("Input matrix contains NaN or Inf values");
  }
  
  arma::mat result(K, T);
  
  // Apply TV prox to each row
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int k = 0; k < K; k++) {
    result.row(k) = prox_tv_condat_1d(X.row(k).t(), lambda_tv).t();
  }
  
  return result;
}

//' Compute Total Variation of Matrix Rows
//' 
//' Computes the sum of total variations across all rows.
//' TV(X) = sum_k sum_t |X[k,t] - X[k,t-1]|
//' 
//' @param X Input matrix (K x T)
//' 
//' @return Total variation value
//' 
//' @export
// [[Rcpp::export]]
double compute_tv_rcpp(const arma::mat& X) {
  // Input validation
  if (X.is_empty()) {
    return 0.0;
  }
  
  if (X.n_cols <= 1) {
    return 0.0;
  }
  
  // Check for numerical issues
  if (X.has_nan() || X.has_inf()) {
    stop("Input matrix contains NaN or Inf values");
  }
  
  return arma::accu(arma::abs(arma::diff(X, 1, 1)));
}

//' Alternative TV Proximal Operator using Dual Method
//' 
//' Implements TV proximal operator using dual formulation.
//' This can be more stable for certain cases.
//' 
//' @param x Input vector
//' @param lambda TV regularization parameter
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' 
//' @return Solution vector
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::vec prox_tv_dual(const arma::vec& x, double lambda,
                       int max_iter = 100, double tol = 1e-8) {
  if (x.is_empty()) {
    stop("Input vector cannot be empty");
  }

  int n = x.n_elem;

  
  if (n <= 2 || lambda <= 0) {


  if (x.has_nan() || x.has_inf()) {
    stop("Input vector contains NaN or Inf values");
  }
  if (!std::isfinite(lambda) || lambda < 0) {
    stop("lambda must be non-negative and finite");
  }
  if (max_iter <= 0) {
    stop("max_iter must be positive");
  }
  if (!std::isfinite(tol) || tol <= 0) {
    stop("tol must be positive and finite");
  }

  if (n <= 1 || lambda == 0) {

    return x;
  }
  
  // Initialize dual variable
  arma::vec p = zeros<vec>(n - 1);
  arma::vec p_old = p;
  
  // Step size for dual updates
  double tau = 0.25;
  
  for (int iter = 0; iter < max_iter; iter++) {
    // Gradient step
    arma::vec grad(n - 1);
    grad(0) = x(0) - x(1) + lambda * (p(0));
    for (int i = 1; i < n - 2; i++) {
      grad(i) = x(i) - x(i+1) + lambda * (p(i) - p(i-1));
    }
    grad(n-2) = x(n-2) - x(n-1) - lambda * p(n-3);
    
    // Update dual variable
    p = p - tau * grad;
    
    // Project onto [-1, 1]
    p = clamp(p, -1.0, 1.0);
    
    // Check convergence
    double change = norm(p - p_old, 2);
    if (change < tol) {
      break;
    }
    
    p_old = p;
  }
  
  // Recover primal solution
  arma::vec z = x;
  z(0) -= lambda * p(0);
  for (int i = 1; i < n - 1; i++) {
    z(i) -= lambda * (p(i) - p(i-1));
  }
  z(n-1) += lambda * p(n-2);
  
  return z;
}
