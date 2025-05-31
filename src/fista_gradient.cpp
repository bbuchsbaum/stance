#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;


// Forward declaration for transposed convolution helper
arma::mat convolve_transpose_rcpp(const arma::mat& X, const arma::vec& hrf,
                                  int n_threads = 0);

arma::mat compute_gradient_fista_precomp_rcpp(const arma::mat& WtY,
                                              const arma::mat& WtW,
                                              const arma::mat& H_star_X,
                                              const arma::vec& hrf_kernel);


//' Compute FISTA Gradient for CLD
//' 
//' Computes the gradient of the least squares term for FISTA optimization.
//' Supports both full computation and efficient pre-computed WtY mode.
//' 
//' @param Y_or_WtY Either Y (V x T) or pre-computed WtY (K x T) matrix
//' @param W Spatial maps matrix (V x K)
//' @param H_star_X Convolved states matrix (K x T)
//' @param hrf_kernel HRF kernel vector
//' @param precomputed_WtY Logical indicating if first argument is WtY
//' @param WtW_precomp Optional pre-computed W'W matrix (K x K)
//' 
//' @return Gradient matrix (K x T)
//' 
//' @export
// [[Rcpp::export]]
arma::mat compute_gradient_fista_rcpp(const arma::mat& Y_or_WtY,
                                      const arma::mat& W,
                                      const arma::mat& H_star_X,
                                      const arma::vec& hrf_kernel,
                                      bool precomputed_WtY = false,
                                      const arma::mat& WtW_precomp = arma::mat()) {
  
  // Input validation
  if (Y_or_WtY.is_empty() || W.is_empty() || H_star_X.is_empty() || hrf_kernel.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }
  
  // Dimension checks
  if (W.n_cols != H_star_X.n_rows) {
    stop("Dimension mismatch: W.n_cols (%d) != H_star_X.n_rows (%d)", 
         W.n_cols, H_star_X.n_rows);
  }
  
  if (precomputed_WtY) {
    if (Y_or_WtY.n_rows != W.n_cols) {
      stop("Dimension mismatch: WtY.n_rows (%d) != W.n_cols (%d)", 
           Y_or_WtY.n_rows, W.n_cols);
    }
    if (Y_or_WtY.n_cols != H_star_X.n_cols) {
      stop("Dimension mismatch: WtY.n_cols (%d) != H_star_X.n_cols (%d)", 
           Y_or_WtY.n_cols, H_star_X.n_cols);
    }
  } else {
    if (Y_or_WtY.n_rows != W.n_rows) {
      stop("Dimension mismatch: Y.n_rows (%d) != W.n_rows (%d)", 
           Y_or_WtY.n_rows, W.n_rows);
    }
    if (Y_or_WtY.n_cols != H_star_X.n_cols) {
      stop("Dimension mismatch: Y.n_cols (%d) != H_star_X.n_cols (%d)", 
           Y_or_WtY.n_cols, H_star_X.n_cols);
    }
  }
  
  arma::mat Grad_term;
  
  if (!precomputed_WtY) {
    // Full computation: Grad_term = W'(Y - W * H_star_X)
    // Residual = Y - W * H_star_X
    arma::mat Residual = Y_or_WtY - W * H_star_X;

    // Grad_term = W' * Residual
    Grad_term = W.t() * Residual;

    // Apply transposed convolution with time-reversed HRF
    arma::mat Grad_L2 = convolve_transpose_rcpp(Grad_term, hrf_kernel);
    return -Grad_L2;

  } else {
    // Efficient computation using pre-computed WtY
    arma::mat WtW;
    if (WtW_precomp.is_empty()) {
      WtW = W.t() * W;
    } else {
      WtW = WtW_precomp;
    }

    // Delegate to the simplified pre-computed interface
    return compute_gradient_fista_precomp_rcpp(Y_or_WtY, WtW,
                                               H_star_X, hrf_kernel);
  }
}

//' Compute FISTA Gradient with Pre-computed Terms
//'
//' Simplified version of \code{compute_gradient_fista_rcpp} that operates
//' directly on pre-computed \eqn{W'Y} and \eqn{W'W} matrices.
//'
//' @param WtY Pre-computed \eqn{W'Y} matrix (K x T)
//' @param WtW Pre-computed \eqn{W'W} matrix (K x K)
//' @param H_star_X Convolved states matrix (K x T)
//' @param hrf_kernel HRF kernel vector
//'
//' @return Gradient matrix (K x T)
//'
//' @export
// [[Rcpp::export]]
arma::mat compute_gradient_fista_precomp_rcpp(const arma::mat& WtY,
                                              const arma::mat& WtW,
                                              const arma::mat& H_star_X,
                                              const arma::vec& hrf_kernel) {
  // Input validation
  if (WtY.is_empty() || WtW.is_empty() || H_star_X.is_empty() ||
      hrf_kernel.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }

  if (WtY.n_rows != WtW.n_rows || WtW.n_rows != WtW.n_cols) {
    stop("Dimension mismatch between WtY and WtW");
  }
  if (WtY.n_rows != H_star_X.n_rows || WtY.n_cols != H_star_X.n_cols) {
    stop("Dimension mismatch between inputs");
  }

  // Gradient term in measurement space
  arma::mat Grad_term = WtY - WtW * H_star_X;

  // Apply transposed convolution with time-reversed HRF

  // This is the gradient with respect to X (before convolution)
  arma::mat Grad_L2 = convolve_transpose_rcpp(Grad_term, hrf_kernel, 0);
  

  return -Grad_L2;
}

//' Transposed Convolution Helper
//' 
//' Performs transposed convolution (correlation) with time-reversed HRF.
//' This is equivalent to convolution with reversed kernel.
//' Supports OpenMP parallelization when available.
//' 
//' @param X Matrix with signals in rows (K x T)
//' @param hrf HRF kernel vector
//' @param n_threads Number of threads to use (0 = auto)
//' 
//' @return Transposed convolution result (K x T)
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat convolve_transpose_rcpp(const arma::mat& X, const arma::vec& hrf,
                                  int n_threads = 0) {
  // Input validation
  if (X.is_empty()) {
    stop("Input matrix X cannot be empty");
  }
  if (hrf.is_empty()) {
    stop("HRF kernel cannot be empty");
  }
  
  int K = X.n_rows;
  int T = X.n_cols;
  int h_len = hrf.n_elem;

#ifdef _OPENMP
  int max_threads = omp_get_max_threads();
  if (n_threads <= 0 || n_threads > max_threads) {
    n_threads = max_threads;
  }
  omp_set_num_threads(n_threads);
#endif
  
  // Check for potential overflow
  if (h_len > T) {
    warning("HRF length (%d) exceeds time series length (%d)", h_len, T);
  }
  
  // Time-reverse the HRF
  arma::vec hrf_rev = reverse(hrf);
  
  // Initialize output
  arma::mat result(K, T, fill::zeros);
  
  // Perform convolution for each row
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (int k = 0; k < K; k++) {
    for (int t = 0; t < T; t++) {
      double sum = 0.0;
      
      // Convolution sum
      for (int tau = 0; tau < h_len; tau++) {
        int t_idx = t - tau;
        if (t_idx >= 0 && t_idx < T) {
          sum += X(k, t_idx) * hrf_rev(tau);
        }
      }
      
      result(k, t) = sum;
    }
  }
  
  return result;
}

//' Compute Lipschitz Constant for FISTA
//' 
//' Estimates the Lipschitz constant of the gradient using power method.
//' L = largest eigenvalue of (H'H ⊗ W'W) where ⊗ is Kronecker product.
//' 
//' @param W Spatial maps matrix (V x K)
//' @param hrf_kernel HRF kernel vector
//' @param max_iter Maximum iterations for power method
//' @param tol Convergence tolerance
//' 
//' @return Estimated Lipschitz constant
//' 
//' @export
// [[Rcpp::export]]
double estimate_lipschitz_rcpp(const arma::mat& W,
                               const arma::vec& hrf_kernel,
                               int max_iter = 30,
                               double tol = 1e-6) {
  
  // Input validation
  if (W.is_empty()) {
    stop("W matrix cannot be empty");
  }
  if (hrf_kernel.is_empty()) {
    stop("HRF kernel cannot be empty");
  }
  
  // Compute W'W
  arma::mat WtW = W.t() * W;
  
  // Check for numerical issues
  if (WtW.has_nan() || WtW.has_inf()) {
    stop("Numerical issues detected in W'W computation");
  }
  
  // For the convolution operator, the spectral norm is bounded by
  // ||H||_2^2 * ||W||_2^2
  // where ||H||_2 is the l2 norm of the HRF
  
  // Compute ||hrf||_2^2
  double hrf_norm_sq = dot(hrf_kernel, hrf_kernel);
  
  // Compute largest eigenvalue of W'W using power method
  int K = WtW.n_rows;
  if (K == 0) {
    stop("W matrix has zero columns");
  }
  
  arma::vec x = randn<vec>(K);
  double x_norm = norm(x, 2);
  if (x_norm < 1e-10) {
    // Extremely unlikely, but handle degenerate case
    x.ones();
    x_norm = norm(x, 2);
  }
  x = x / x_norm;
  
  double lambda_old = 0.0;
  double lambda_new = 0.0;
  
  for (int iter = 0; iter < max_iter; iter++) {
    // Apply operator
    arma::vec Ax = WtW * x;
    
    // Compute Rayleigh quotient
    lambda_new = dot(x, Ax);
    
    // Check convergence
    if (std::abs(lambda_new - lambda_old) < tol) {
      break;
    }
    
    // Update
    double Ax_norm = norm(Ax, 2);
    if (Ax_norm < 1e-10) {
      // Matrix is near-singular
      warning("W'W appears to be near-singular");
      break;
    }
    x = Ax / Ax_norm;
    lambda_old = lambda_new;
  }
  
  // Return L = ||H||_2^2 * lambda_max(W'W)
  return hrf_norm_sq * lambda_new;
}

