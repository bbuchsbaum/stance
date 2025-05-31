#' CLD Mathematical Helper Functions
#'
#' Internal mathematical utilities for the Continuous Linear Decoder.
#' These functions provide efficient implementations of key operations
#' used throughout the CLD algorithm.
#'
#' @name cld_math
#' @keywords internal
NULL

#' Fast Matrix Multiplication with Caching
#'
#' Computes matrix products with optional caching of intermediate results.
#'
#' @param A First matrix
#' @param B Second matrix
#' @param cache_key Optional key for caching result
#' @param cache_env Environment for storing cached results
#' 
#' @return Matrix product A %*% B
#' @keywords internal
fast_matmul <- function(A, B, cache_key = NULL, cache_env = NULL) {
  if (!is.null(cache_key) && !is.null(cache_env)) {
    if (exists(cache_key, envir = cache_env)) {
      return(get(cache_key, envir = cache_env))
    }
  }
  
  result <- A %*% B
  
  if (!is.null(cache_key) && !is.null(cache_env)) {
    assign(cache_key, result, envir = cache_env)
  }
  
  result
}

#' Compute Row-wise Convolution
#'
#' Helper function for efficient row-wise convolution operations.
#'
#' @param X Matrix with signals in rows
#' @param h Kernel vector
#' @param method Convolution method ("direct" or "fft")
#' 
#' @return Convolved matrix
#' @keywords internal
convolve_rows <- function(X, h, method = "direct") {
  convolve_with_hrf(X, h, method = method)
}

#' Soft Thresholding Operator
#'
#' Applies soft thresholding for L1 regularization.
#'
#' @param x Input vector or matrix
#' @param lambda Threshold parameter
#' 
#' @return Soft-thresholded result
#' @keywords internal
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

#' Compute Total Variation
#'
#' Computes the total variation (L1 norm of differences) along rows.
#'
#' @param X Matrix (signals in rows)
#' 
#' @return Total variation value
#' @keywords internal
compute_tv <- function(X) {
  if (ncol(X) <= 1) return(0)
  
  # Compute differences
  diffs <- X[, -1] - X[, -ncol(X)]
  
  # L1 norm
  sum(abs(diffs))
}

#' Compute Objective Function Value
#'
#' Evaluates the CLD objective function.
#'
#' @param Y Data matrix (V x T)
#' @param W Spatial maps (V x K)
#' @param X States (K x T)
#' @param hrf HRF kernel
#' @param lambda_tv TV regularization parameter
#' 
#' @return Objective function value
#' @keywords internal
compute_objective <- function(Y, W, X, hrf, lambda_tv) {
  # Convolve states with HRF
  HX <- convolve_rows(X, hrf)
  
  # Reconstruction error
  residual <- Y - W %*% HX
  reconstruction_error <- 0.5 * sum(residual^2)
  
  # TV penalty
  tv_penalty <- lambda_tv * compute_tv(X)
  
  reconstruction_error + tv_penalty
}

#' Project onto Orthogonal Columns
#'
#' Projects a matrix to have orthonormal columns.
#'
#' @param U Matrix to project
#' 
#' @return Matrix with orthonormal columns
#' @keywords internal
project_orthogonal <- function(U) {
  svd_result <- svd(U, nu = ncol(U), nv = 0)
  svd_result$u
}

#' Block-wise Processing
#'
#' Processes large matrices in blocks for memory efficiency.
#'
#' @param X Large matrix
#' @param fun Function to apply to each block
#' @param block_size Number of rows per block
#' @param ... Additional arguments to fun
#' 
#' @return Processed matrix
#' @keywords internal
process_blocks <- function(X, fun, block_size = 1000, ...) {
  n_rows <- nrow(X)
  n_blocks <- ceiling(n_rows / block_size)
  
  results <- list()
  
  for (b in seq_len(n_blocks)) {
    start_idx <- (b - 1) * block_size + 1
    end_idx <- min(b * block_size, n_rows)
    
    block <- X[start_idx:end_idx, , drop = FALSE]
    results[[b]] <- fun(block, ...)
  }
  
  do.call(rbind, results)
}

#' Estimate Operator Norm via Power Method
#'
#' Estimates the largest singular value of a linear operator.
#'
#' @param A Linear operator (function or matrix)
#' @param dim Dimension for random initialization
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance
#' 
#' @return Estimated operator norm
#' @keywords internal
estimate_operator_norm <- function(A, dim, max_iter = 30, tol = 1e-6) {
  # Initialize with random vector
  x <- rnorm(dim)
  x <- x / sqrt(sum(x^2))
  
  lambda_old <- 0
  
  for (iter in seq_len(max_iter)) {
    # Apply operator
    if (is.function(A)) {
      Ax <- A(x)
    } else {
      Ax <- A %*% x
    }
    
    # Compute Rayleigh quotient
    lambda <- sqrt(sum(Ax^2))
    
    # Check convergence
    if (abs(lambda - lambda_old) < tol) {
      break
    }
    
    # Update
    x <- Ax / lambda
    lambda_old <- lambda
  }
  
  lambda
}

#' Initialize State Matrix
#'
#' Creates initial state matrix for optimization.
#'
#' @param K Number of states
#' @param T Number of time points
#' @param method Initialization method ("zeros", "random", "uniform")
#' 
#' @return K x T initial state matrix
#' @keywords internal
initialize_states <- function(K, T, method = "zeros") {
  switch(method,
    "zeros" = matrix(0, K, T),
    "random" = matrix(rnorm(K * T, sd = 0.1), K, T),
    "uniform" = matrix(1/K, K, T),
    stop("Unknown initialization method: ", method)
  )
}

#' Check Convergence
#'
#' Checks if optimization has converged based on relative change.
#'
#' @param values Vector of objective values
#' @param tol Relative tolerance
#' @param window Number of iterations to check
#' 
#' @return Logical indicating convergence
#' @keywords internal
check_convergence <- function(values, tol = 1e-4, window = 5) {
  n <- length(values)
  
  if (n < window) {
    return(FALSE)
  }
  
  recent <- values[(n - window + 1):n]
  rel_change <- abs(recent[window] - recent[1]) / (abs(recent[1]) + 1e-10)
  
  rel_change < tol
}

#' Create Progress Bar
#'
#' Creates a simple text progress bar for iterations.
#'
#' @param total Total number of iterations
#' @param width Width of progress bar
#' 
#' @return Function to update progress
#' @keywords internal
create_progress_bar <- function(total, width = 50) {
  function(current) {
    if (current == 1) {
      cat("Progress: |", rep(" ", width), "| 0%\r", sep = "")
    }
    
    pct <- current / total
    filled <- floor(pct * width)
    
    cat("Progress: |", 
        rep("=", filled), 
        rep(" ", width - filled), 
        "| ", 
        sprintf("%.0f%%", pct * 100),
        "\r", 
        sep = "")
    
    if (current == total) {
      cat("\n")
    }
    
    flush.console()
  }
}
