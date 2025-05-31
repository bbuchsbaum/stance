#' HRF Utilities for stance
#'
#' Functions for HRF processing leveraging fmrireg infrastructure.
#' Provides efficient convolution operations and basis function handling
#' for both fixed HRFs (CLD) and voxel-wise HRFs (CBD).
#'
#' @name hrf_utils
NULL

#' Setup HRF Kernel
#'
#' Accepts fmrireg HRF specifications or numeric vectors and normalizes
#' them to unit energy for consistent convolution operations.
#'
#' @param spec HRF specification: can be a character string matching fmrireg
#'   HRF types (e.g., "spmg1", "spmg2", "bspline"), a numeric vector, or
#'   an fmrireg HRF object
#' @param TR Repetition time in seconds (default: 2)
#' @param len Length of HRF in seconds (default: 32)
#' @param normalize Logical, whether to normalize to unit energy (default: TRUE)
#' 
#' @return Numeric vector representing the HRF kernel
#' 
#' @export
#' @examples
#' \dontrun{
#' # Using fmrireg canonical HRF
#' hrf1 <- setup_hrf_kernel("spmg1", TR = 2)
#' 
#' # Using custom numeric vector
#' custom_hrf <- c(0, 0.1, 0.5, 0.8, 0.6, 0.3, 0.1, 0)
#' hrf2 <- setup_hrf_kernel(custom_hrf)
#' }
setup_hrf_kernel <- function(spec = "spmg1", TR = 2, len = 32, normalize = TRUE) {
  # Handle numeric vector input
  if (is.numeric(spec)) {
    hrf_vec <- spec
    if (normalize && sum(abs(hrf_vec)) > 0) {
      # Normalize to unit energy (sum of squares = 1)
      hrf_vec <- hrf_vec / sqrt(sum(hrf_vec^2))
    }
    return(hrf_vec)
  }
  
  # Handle fmrireg HRF specifications
  if (is.character(spec)) {
    # Map common names to fmrireg constants
    spec_map <- c(
      "canonical" = "spmg1",
      "spm" = "spmg1",
      "spmg1" = "spmg1",
      "spmg2" = "spmg2",
      "spmg3" = "spmg3",
      "gamma" = "gamma",
      "gaussian" = "gaussian",
      "bspline" = "bspline",
      "fir" = "fir"
    )
    
    spec_lower <- tolower(spec)
    if (spec_lower %in% names(spec_map)) {
      spec <- spec_map[spec_lower]
    }
    
    # Generate HRF using fmrireg
    # Create sampling frame for HRF evaluation
    sframe <- fmrireg::sampling_frame(1, len, TR)
    
    # Get appropriate HRF basis from fmrireg
    hrf_basis <- switch(spec,
      "spmg1" = fmrireg::HRF_SPMG1,
      "spmg2" = fmrireg::HRF_SPMG2,
      "spmg3" = fmrireg::HRF_SPMG3,
      "gamma" = fmrireg::HRF_GAMMA,
      "gaussian" = fmrireg::HRF_GAUSSIAN,
      "bspline" = fmrireg::HRF_BSPLINE,
      "fir" = fmrireg::HRF_FIR,
      stop("Unknown HRF specification: ", spec)
    )
    
    # Evaluate HRF at time points
    hrf_vec <- fmrireg::evaluate(hrf_basis, seq(0, len, by = TR))
    
    # For multi-column bases, take the first (canonical) column
    if (is.matrix(hrf_vec)) {
      hrf_vec <- hrf_vec[, 1]
    }
    
    if (normalize && sum(abs(hrf_vec)) > 0) {
      hrf_vec <- hrf_vec / sqrt(sum(hrf_vec^2))
    }
    
    return(hrf_vec)
  }
  
  # Handle fmrireg HRF objects directly
  if (inherits(spec, "hrf")) {
    hrf_vec <- fmrireg::evaluate(spec, seq(0, len, by = TR))
    if (is.matrix(hrf_vec)) {
      hrf_vec <- hrf_vec[, 1]
    }
    if (normalize && sum(abs(hrf_vec)) > 0) {
      hrf_vec <- hrf_vec / sqrt(sum(hrf_vec^2))
    }
    return(hrf_vec)
  }
  
  stop("HRF specification must be a character string, numeric vector, or fmrireg HRF object")
}

#' Convolve with HRF
#'
#' Efficient convolution with automatic FFT selection based on signal length.
#' Handles both single HRF (CLD) and multiple HRFs (CBD) scenarios.
#'
#' @param X Matrix with signals in rows (K x T)
#' @param hrf HRF kernel vector or matrix (for multiple HRFs)
#' @param use_fft_threshold Integer, use FFT when T > threshold (default: 256)
#' @param method Character, convolution method ("auto", "direct", "fft")
#' 
#' @return Matrix of convolved signals (same dimensions as X)
#' 
#' @export
convolve_with_hrf <- function(X, hrf, use_fft_threshold = 256, method = "auto") {
  # Validate inputs
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  n_signals <- nrow(X)
  n_time <- ncol(X)
  
  # Determine convolution method
  if (method == "auto") {
    method <- ifelse(n_time > use_fft_threshold, "fft", "direct")
  }
  
  # Single HRF case
  if (is.vector(hrf)) {
    if (method == "fft") {
      # FFT-based convolution
      result <- t(apply(X, 1, function(row) {
        convolve_fft(row, hrf)
      }))
    } else {
      # Direct convolution
      result <- t(apply(X, 1, function(row) {
        conv_out <- stats::convolve(row, rev(hrf), type = "open")
        # Trim to original length
        conv_out[seq_len(n_time)]
      }))
    }
  } else if (is.matrix(hrf)) {
    # Multiple HRFs case (for future CBD support)
    if (nrow(hrf) != n_signals) {
      stop("Number of HRFs must match number of signals")
    }
    
    result <- matrix(0, n_signals, n_time)
    for (i in seq_len(n_signals)) {
      if (method == "fft") {
        result[i, ] <- convolve_fft(X[i, ], hrf[i, ])
      } else {
        conv_out <- stats::convolve(X[i, ], rev(hrf[i, ]), type = "open")
        result[i, ] <- conv_out[seq_len(n_time)]
      }
    }
  } else {
    stop("HRF must be a vector or matrix")
  }
  
  result
}

#' FFT-based Convolution
#'
#' Internal function for FFT-based convolution.
#'
#' @param x Signal vector
#' @param h Kernel vector
#' @return Convolved signal (trimmed to original length)
#' @keywords internal
convolve_fft <- function(x, h) {
  n_x <- length(x)
  n_h <- length(h)
  n_fft <- 2^ceiling(log2(n_x + n_h - 1))
  
  # Zero-pad signals
  x_pad <- c(x, rep(0, n_fft - n_x))
  h_pad <- c(h, rep(0, n_fft - n_h))
  
  # FFT convolution
  X_fft <- stats::fft(x_pad)
  H_fft <- stats::fft(h_pad)
  Y_fft <- X_fft * H_fft
  
  # Inverse FFT and trim
  y <- Re(stats::fft(Y_fft, inverse = TRUE)) / n_fft
  y[seq_len(n_x)]
}

#' Generate HRF Basis Matrix
#'
#' Creates basis matrices for flexible HRF modeling, supporting both
#' fixed basis sets (CLD) and voxel-wise expansions (CBD).
#'
#' @param basis_type Character string specifying basis type
#' @param TR Repetition time in seconds
#' @param len Length of HRF in seconds
#' @param n_basis Number of basis functions (for FIR, B-spline)
#' 
#' @return Matrix with basis functions in columns
#' 
#' @export
hrf_basis_matrix <- function(basis_type = "spmg3", TR = 2, len = 32, n_basis = NULL) {
  # Create sampling frame
  sframe <- fmrireg::sampling_frame(1, len, TR)
  time_points <- seq(0, len, by = TR)
  
  # Get appropriate basis
  basis_obj <- switch(tolower(basis_type),
    "spmg1" = fmrireg::HRF_SPMG1,
    "spmg2" = fmrireg::HRF_SPMG2,
    "spmg3" = fmrireg::HRF_SPMG3,
    "gamma" = fmrireg::HRF_GAMMA,
    "gaussian" = fmrireg::HRF_GAUSSIAN,
    "bspline" = {
      if (is.null(n_basis)) n_basis <- 5
      fmrireg::HRF_BSPLINE
    },
    "fir" = {
      if (is.null(n_basis)) n_basis <- 12
      fmrireg::HRF_FIR
    },
    stop("Unknown basis type: ", basis_type)
  )
  
  # Evaluate basis
  basis_mat <- fmrireg::evaluate(basis_obj, time_points)
  
  # Ensure matrix output
  if (!is.matrix(basis_mat)) {
    basis_mat <- matrix(basis_mat, ncol = 1)
  }
  
  # Add attributes
  attr(basis_mat, "basis_type") <- basis_type
  attr(basis_mat, "TR") <- TR
  attr(basis_mat, "time_points") <- time_points
  
  basis_mat
}

#' Create HRF Basis (neuroim2 integration)
#'
#' Convenience wrapper that generates an HRF basis matrix using
#' `hrf_basis_matrix()` and returns it for use with neuroim2 based
#' workflows.
#'
#' @param type Character string specifying the basis type
#' @param TR Repetition time in seconds
#' @param duration Total duration of the HRF basis in seconds
#'
#' @return Matrix of basis functions (time points x basis functions)
#' @export
create_hrf_basis_neuroim2 <- function(type = "spmg3", TR = 2, duration = 32) {
  validate_hrf_spec(type, context = "hrf basis type")
  hrf_basis_matrix(basis_type = type, TR = TR, len = duration)
}

#' Validate HRF Specification
#'
#' Checks HRF specifications and provides informative error messages.
#'
#' @param hrf_spec HRF specification to validate
#' @param context Character string describing the context (for error messages)
#' 
#' @return TRUE if valid, otherwise stops with error
#' 
#' @export
validate_hrf_spec <- function(hrf_spec, context = "HRF specification") {
  # Check for NULL
  if (is.null(hrf_spec)) {
    stop(context, " cannot be NULL")
  }
  
  # Valid character specifications
  valid_specs <- c("spmg1", "spmg2", "spmg3", "gamma", "gaussian", 
                   "bspline", "fir", "canonical", "spm")
  
  if (is.character(hrf_spec)) {
    if (!tolower(hrf_spec) %in% valid_specs) {
      stop(context, " '", hrf_spec, "' is not a recognized HRF type. ",
           "Valid options: ", paste(valid_specs, collapse = ", "))
    }
  } else if (is.numeric(hrf_spec)) {
    if (length(hrf_spec) < 2) {
      stop(context, " as numeric vector must have length >= 2")
    }
    if (any(!is.finite(hrf_spec))) {
      stop(context, " contains non-finite values")
    }
    if (sum(abs(hrf_spec)) == 0) {
      stop(context, " is all zeros")
    }
  } else if (!inherits(hrf_spec, "hrf")) {
    stop(context, " must be a character string, numeric vector, or fmrireg HRF object")
  }
  
  TRUE
}

#' Create HRF Convolution Matrix
#'
#' Creates a sparse Toeplitz matrix for HRF convolution operations.
#' Useful for CBD implementations where explicit matrix operations are needed.
#'
#' @param hrf HRF kernel vector
#' @param n_time Number of time points
#' @param sparse Logical, whether to return sparse matrix (default: TRUE)
#' 
#' @return Convolution matrix (n_time x n_time)
#' 
#' @export
hrf_convolution_matrix <- function(hrf, n_time, sparse = TRUE) {
  n_hrf <- length(hrf)
  
  if (sparse) {
    # Create sparse Toeplitz matrix
    # First column: [hrf, zeros]
    # First row: [hrf[1], zeros]
    i <- rep(1:n_time, each = n_hrf)
    j <- unlist(lapply(1:n_time, function(t) {
      idx <- (t - n_hrf + 1):t
      idx[idx > 0 & idx <= n_time]
    }))
    x <- rep(rev(hrf), n_time)[seq_along(j)]
    
    # Remove entries where j would be out of bounds
    valid <- j > 0 & j <= n_time
    
    Matrix::sparseMatrix(
      i = i[valid],
      j = j[valid],
      x = x[valid],
      dims = c(n_time, n_time)
    )
  } else {
    # Dense Toeplitz matrix
    conv_mat <- matrix(0, n_time, n_time)
    for (t in 1:n_time) {
      for (tau in 1:min(n_hrf, t)) {
        if (t - tau + 1 > 0) {
          conv_mat[t, t - tau + 1] <- hrf[tau]
        }
      }
    }
    conv_mat
  }
}

#' Transposed Convolution with HRF
#'
#' Performs transposed convolution (correlation) with time-reversed HRF.
#' Used in gradient computations for FISTA optimization.
#'
#' @param X Matrix with signals in rows (K x T)
#' @param hrf HRF kernel vector
#' @param method Character, convolution method ("auto", "direct", "fft")
#' 
#' @return Matrix of transposed-convolved signals
#' 
#' @export
convolve_with_hrf_transposed <- function(X, hrf, method = "auto") {
  # Time-reverse the HRF for transposed convolution
  hrf_reversed <- rev(hrf)
  
  # Use the same convolution function with reversed kernel
  convolve_with_hrf(X, hrf_reversed, method = method)
}
