#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

//' Parallel Voxel-wise GLM Fitting
//' 
//' Fits GLMs for multiple voxels in parallel using OpenMP.
//' Each voxel's time series is regressed on the convolved design matrix.
//' 
//' @param Y Data matrix (V x T)
//' @param X_conv Convolved design matrix (T x K)
//' @param n_threads Number of threads to use (0 = auto)
//' @param chunk_size Size of chunks for dynamic scheduling
//' 
//' @return Beta coefficients matrix (V x K)
//' 
//' @export
// [[Rcpp::export]]
arma::mat parallel_glm_fit_rcpp(const arma::mat& Y,
                                const arma::mat& X_conv,
                                int n_threads = 0,
                                int chunk_size = 100) {
  
  // Input validation
  if (Y.is_empty() || X_conv.is_empty()) {
    stop("Input matrices cannot be empty");
  }
  
  int V = Y.n_rows;
  int T = Y.n_cols;
  int K = X_conv.n_cols;
  
  // Dimension checks
  if (X_conv.n_rows != T) {
    stop("Dimension mismatch: X_conv.n_rows (%d) != Y.n_cols (%d)", 
         X_conv.n_rows, T);
  }
  
  // Initialize output
  arma::mat B(V, K, fill::zeros);
  
  // Set number of threads
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n_threads <= 0 || n_threads > max_threads) {
    n_threads = max_threads;
  }
  omp_set_num_threads(n_threads);
  
  if (V > 1000) {  // Only report for large problems
    Rcout << "Using " << n_threads << " threads for GLM fitting\n";
  }
#endif
  
  // Pre-compute X'X and its inverse for efficiency
  arma::mat XtX = X_conv.t() * X_conv;
  
  // Check for singularity
  double det_XtX = det(XtX);
  if (std::abs(det_XtX) < 1e-10) {
    stop("Design matrix is singular or near-singular");
  }
  
  arma::mat XtX_inv;
  bool inv_success = inv_sympd(XtX_inv, XtX);
  if (!inv_success) {
    // Fallback to pseudo-inverse
    warning("Using pseudo-inverse due to numerical issues in X'X");
    XtX_inv = pinv(XtX);
  }
  
  // Pre-compute (X'X)^(-1) * X' for efficiency
  arma::mat XtX_inv_Xt = XtX_inv * X_conv.t();
  
  // Parallel loop over voxels
  #pragma omp parallel for schedule(dynamic, chunk_size)
  for (int v = 0; v < V; v++) {
    // Extract voxel time series
    arma::vec y_v = Y.row(v).t();
    
    // Check for missing data
    if (y_v.has_nan() || y_v.has_inf()) {
      // Skip this voxel, leave coefficients as zeros
      continue;
    }
    
    // Compute beta = (X'X)^(-1) * X' * y
    arma::vec beta_v = XtX_inv_Xt * y_v;
    
    // Store in output matrix
    B.row(v) = beta_v.t();
  }
  
  return B;
}

//' Parallel GLM with Regularization
//' 
//' Fits regularized GLMs (Ridge regression) for multiple voxels in parallel.
//' 
//' @param Y Data matrix (V x T)
//' @param X_conv Convolved design matrix (T x K)
//' @param lambda Ridge penalty parameter
//' @param n_threads Number of threads to use (0 = auto)
//' 
//' @return Beta coefficients matrix (V x K)
//' 
//' @export
// [[Rcpp::export]]
arma::mat parallel_ridge_glm_rcpp(const arma::mat& Y,
                                  const arma::mat& X_conv,
                                  double lambda = 0.01,
                                  int n_threads = 0) {
  
  // Input validation
  if (Y.is_empty() || X_conv.is_empty()) {
    stop("Input matrices cannot be empty");
  }
  if (!std::isfinite(lambda) || lambda < 0) {
    stop("lambda must be non-negative and finite");
  }
  
  int V = Y.n_rows;
  int T = Y.n_cols;
  int K = X_conv.n_cols;
  
  // Dimension checks
  if (X_conv.n_rows != T) {
    stop("Dimension mismatch: X_conv.n_rows (%d) != Y.n_cols (%d)", 
         X_conv.n_rows, T);
  }
  
  // Initialize output
  arma::mat B(V, K, fill::zeros);
  
  // Set number of threads
#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n_threads <= 0 || n_threads > max_threads) {
    n_threads = max_threads;
  }
  omp_set_num_threads(n_threads);
#endif
  
  // Pre-compute X'X + lambda*I and its inverse
  arma::mat XtX = X_conv.t() * X_conv;
  arma::mat XtX_reg = XtX + lambda * eye<mat>(K, K);
  
  arma::mat XtX_reg_inv;
  bool inv_success = inv_sympd(XtX_reg_inv, XtX_reg);
  if (!inv_success) {
    stop("Failed to invert regularized X'X matrix");
  }
  
  // Pre-compute (X'X + lambda*I)^(-1) * X'
  arma::mat XtX_reg_inv_Xt = XtX_reg_inv * X_conv.t();
  
  // Parallel loop over voxels
  #pragma omp parallel for schedule(dynamic, 100)
  for (int v = 0; v < V; v++) {
    // Extract voxel time series
    arma::vec y_v = Y.row(v).t();
    
    // Check for missing data
    if (y_v.has_nan() || y_v.has_inf()) {
      continue;
    }
    
    // Compute beta = (X'X + lambda*I)^(-1) * X' * y
    arma::vec beta_v = XtX_reg_inv_Xt * y_v;
    
    // Store in output matrix
    B.row(v) = beta_v.t();
  }
  
  return B;
}

//' Check OpenMP Support
//' 
//' Returns information about OpenMP availability and configuration.
//' 
//' @return List with OpenMP information
//' 
//' @export
// [[Rcpp::export]]
List check_openmp_support() {
  bool has_openmp = false;
  int max_threads = 1;
  int num_procs = 1;
  
#ifdef _OPENMP
  has_openmp = true;
  max_threads = omp_get_max_threads();
  num_procs = omp_get_num_procs();
#endif
  
  return List::create(
    Named("available") = has_openmp,
    Named("max_threads") = max_threads,
    Named("num_processors") = num_procs
  );
}
