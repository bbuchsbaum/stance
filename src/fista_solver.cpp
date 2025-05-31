#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Forward declarations of functions from other files
arma::mat compute_gradient_fista_rcpp(const arma::mat& Y_or_WtY,
                                      const arma::mat& W,
                                      const arma::mat& H_star_X,
                                      const arma::vec& hrf_kernel,
                                      bool precomputed_WtY,
                                      const arma::mat& WtW_precomp = arma::mat());

arma::mat compute_gradient_fista_precomp_rcpp(const arma::mat& WtY,
                                              const arma::mat& WtW,
                                              const arma::mat& H_star_X,
                                              const arma::vec& hrf_kernel);

arma::mat prox_tv_condat_rcpp(const arma::mat& X, double lambda_tv);

arma::mat convolve_rows_rcpp(const arma::mat& X, const arma::vec& hrf);

double compute_tv_rcpp(const arma::mat& X);

//' Fast Row-wise Convolution
//' 
//' Performs convolution of each row with the HRF kernel.
//' 
//' @param X Matrix with signals in rows (K x T)
//' @param hrf HRF kernel vector
//' 
//' @return Convolved matrix (K x T)
//' 
//' @keywords internal
// [[Rcpp::export]]
arma::mat convolve_rows_rcpp(const arma::mat& X, const arma::vec& hrf) {
  // Input validation
  if (X.is_empty() || hrf.is_empty()) {
    stop("Input matrix or HRF kernel cannot be empty");
  }
  
  int K = X.n_rows;
  int T = X.n_cols;
  int h_len = hrf.n_elem;
  
  // Check for potential issues
  if (h_len > T) {
    warning("HRF length (%d) exceeds time series length (%d)", h_len, T);
  }
  
  arma::mat result(K, T, fill::zeros);
  
  // Perform convolution for each row
  for (int k = 0; k < K; k++) {
    for (int t = 0; t < T; t++) {
      double sum = 0.0;
      
      // Convolution sum
      for (int tau = 0; tau < std::min(h_len, t + 1); tau++) {
        sum += X(k, t - tau) * hrf(tau);
      }
      
      result(k, t) = sum;
    }
  }
  
  return result;
}

//' FISTA Solver for CLD
//' 
//' Implements the Fast Iterative Shrinkage-Thresholding Algorithm (FISTA)
//' with Total Variation regularization for continuous state estimation.
//' 
//' @param WtY Pre-computed W'Y matrix (K x T)
//' @param W Spatial maps matrix (V x K)
//' @param hrf_kernel HRF kernel vector
//' @param lambda_tv TV regularization parameter
//' @param L_fista Lipschitz constant
//' @param X_init Initial state matrix (K x T)
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param verbose Print progress
//' 
//' @return List containing:
//'   - X_hat: Estimated states (K x T)
//'   - objective_values: Vector of objective values
//'   - converged: Logical indicating convergence
//'   - iterations: Number of iterations performed
//' 
//' @export
// [[Rcpp::export]]
List fista_tv_rcpp(const arma::mat& WtY,
                   const arma::mat& W,
                   const arma::vec& hrf_kernel,
                   double lambda_tv,
                   double L_fista,
                   const arma::mat& X_init,
                   int max_iter = 100,
                   double tol = 1e-4,
                   bool verbose = false) {
  
  // Input validation
  if (WtY.is_empty() || W.is_empty() || hrf_kernel.is_empty() || X_init.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }
  
  // Dimension checks
  if (WtY.n_rows != W.n_cols) {
    stop("Dimension mismatch: WtY.n_rows (%d) != W.n_cols (%d)", 
         WtY.n_rows, W.n_cols);
  }
  if (WtY.n_rows != X_init.n_rows || WtY.n_cols != X_init.n_cols) {
    stop("Dimension mismatch between WtY and X_init");
  }
  
  // Parameter validation
  if (!std::isfinite(lambda_tv) || lambda_tv < 0) {
    stop("lambda_tv must be non-negative and finite");
  }
  if (!std::isfinite(L_fista) || L_fista <= 0) {
    stop("L_fista must be positive and finite");
  }
  if (max_iter <= 0) {
    stop("max_iter must be positive");
  }
  if (!std::isfinite(tol) || tol <= 0) {
    stop("tol must be positive and finite");
  }
  
  int K = X_init.n_rows;
  int T = X_init.n_cols;
  
  // Initialize variables
  arma::mat X = X_init;
  arma::mat X_old = X;
  arma::mat Z = X;
  double t = 1.0;
  double t_old = 1.0;
  
  // Step size
  double step_size = 1.0 / L_fista;
  if (!std::isfinite(step_size)) {
    stop("Invalid step size computed from L_fista");
  }
  
  // Track objective values
  std::vector<double> objective_values;
  objective_values.reserve(max_iter);
  
  // Pre-compute W'W for efficiency
  arma::mat WtW = W.t() * W;
  
  bool converged = false;
  int iter;
  
  if (verbose) {
    Rcout << "Starting FISTA optimization...\n";
    Rcout << "Step size: " << step_size << "\n";
  }
  
  for (iter = 0; iter < max_iter; iter++) {
    // Store old values
    X_old = X;
    t_old = t;
    
    // Compute gradient at Z
    arma::mat H_star_Z = convolve_rows_rcpp(Z, hrf_kernel);
    arma::mat gradient = compute_gradient_fista_precomp_rcpp(WtY, WtW,
                                                             H_star_Z,
                                                             hrf_kernel);
    
    // Gradient step
    arma::mat X_tilde = Z - step_size * gradient;
    
    // Proximal step (TV denoising)
    X = prox_tv_condat_rcpp(X_tilde, step_size * lambda_tv);
    
    // Update momentum parameter
    t = (1.0 + std::sqrt(1.0 + 4.0 * t_old * t_old)) / 2.0;
    
    // Extrapolation step
    double momentum = (t_old - 1.0) / t;
    Z = X + momentum * (X - X_old);
    
    // Compute objective value (optional, for convergence check)
    if (iter % 10 == 0 || iter < 10) {
      // Compute residual
      arma::mat H_star_X = convolve_rows_rcpp(X, hrf_kernel);
      arma::mat residual = WtY - WtW * H_star_X;
      double reconstruction_error = 0.5 * accu(square(residual));
      
      // Add TV penalty
      double tv_penalty = lambda_tv * compute_tv_rcpp(X);
      
      double obj_val = reconstruction_error + tv_penalty;
      objective_values.push_back(obj_val);
      
      if (verbose && iter % 10 == 0) {
        Rcout << "Iteration " << iter << ": objective = " << obj_val << "\n";
      }
      
      // Check convergence
      if (objective_values.size() > 5) {
        int n = objective_values.size();
        double recent_change = std::abs(objective_values[n-1] - objective_values[n-5]) / 
                               (std::abs(objective_values[n-5]) + 1e-10);
        
        if (recent_change < tol) {
          converged = true;
          if (verbose) {
            Rcout << "Converged at iteration " << iter << "\n";
          }
          break;
        }
      }
    }
  }
  
  if (verbose && !converged) {
    Rcout << "FISTA did not converge within " << max_iter << " iterations\n";
  }
  
  return List::create(
    Named("X_hat") = X,
    Named("objective_values") = objective_values,
    Named("converged") = converged,
    Named("iterations") = iter + 1
  );
}

//' Compute CLD Objective Function
//' 
//' Evaluates the complete objective function including reconstruction error
//' and TV penalty.
//' 
//' @param Y Data matrix (V x T)
//' @param W Spatial maps (V x K)
//' @param X States (K x T)
//' @param hrf HRF kernel
//' @param lambda_tv TV regularization parameter
//' 
//' @return Objective function value
//' 
//' @export
// [[Rcpp::export]]
double compute_objective_rcpp(const arma::mat& Y,
                             const arma::mat& W,
                             const arma::mat& X,
                             const arma::vec& hrf,
                             double lambda_tv) {
  // Input validation
  if (Y.is_empty() || W.is_empty() || X.is_empty() || hrf.is_empty()) {
    stop("Input matrices/vectors cannot be empty");
  }
  
  // Dimension checks
  if (Y.n_rows != W.n_rows) {
    stop("Dimension mismatch: Y.n_rows (%d) != W.n_rows (%d)", 
         Y.n_rows, W.n_rows);
  }
  if (W.n_cols != X.n_rows) {
    stop("Dimension mismatch: W.n_cols (%d) != X.n_rows (%d)", 
         W.n_cols, X.n_rows);
  }
  if (Y.n_cols != X.n_cols) {
    stop("Dimension mismatch: Y.n_cols (%d) != X.n_cols (%d)", 
         Y.n_cols, X.n_cols);
  }
  
  // Convolve states with HRF
  arma::mat HX = convolve_rows_rcpp(X, hrf);
  
  // Reconstruction error
  arma::mat residual = Y - W * HX;
  double reconstruction_error = 0.5 * accu(square(residual));
  
  // TV penalty
  double tv_penalty = lambda_tv * compute_tv_rcpp(X);
  
  return reconstruction_error + tv_penalty;
}
