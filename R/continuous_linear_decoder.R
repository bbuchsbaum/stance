#' Continuous Linear Decoder (CLD)
#'
#' R6 class implementing the Continuous Linear Decoder for fMRI analysis.
#' Provides fast, deterministic continuous state estimation using GLM+SVD
#' for spatial map learning and FISTA with Total Variation regularization
#' for temporal state estimation.
#'
#' @description
#' The CLD solves: \deqn{\hat{X} = \arg\min_{X} ||Y - W(H * X)||_F^2 + \lambda_{TV} ||∇_t X||_1}
#' where Y is the fMRI data, W are spatial maps, H is the HRF, X are state activations,
#' and the Total Variation penalty enforces temporal smoothness.
#'
#' @export
#' @importFrom R6 R6Class
#' @importFrom Matrix sparseMatrix
#' @importFrom stats lm convolve
ContinuousLinearDecoder <- R6::R6Class(
  "ContinuousLinearDecoder",
  
  public = list(
    #' @description
    #' Create a new ContinuousLinearDecoder object.
    #'
    #' @param Y Input fMRI data (V x T matrix or NeuroVec object)
    #' @param S_design Design matrix for supervised learning (K x T)
    #' @param hrf HRF specification (see setup_hrf_kernel)
    #' @param rank Rank for low-rank approximation
    #' @param lambda_tv Total variation regularization parameter
    #' @param verbose Logical, print progress messages
    #'
    #' @return A new ContinuousLinearDecoder object
    initialize = function(Y, S_design, hrf = "spmg1", rank = 20, 
                          lambda_tv = 0.01, verbose = FALSE) {
      # Validate and process input data
      if (verbose) cat("Initializing Continuous Linear Decoder...\n")
      
      # Validate fMRI input
      private$.Y_info <- validate_fmri_input(Y, verbose = verbose)
      private$.Y <- private$.Y_info$data
      V <- private$.Y_info$dims["V"]
      T <- private$.Y_info$dims["T"]
      
      # Validate S_design
      if (!is.matrix(S_design)) {
        stop("S_design must be a matrix")
      }
      if (ncol(S_design) != T) {
        stop(sprintf("S_design has %d columns but Y has %d timepoints", 
                     ncol(S_design), T))
      }
      private$.S_design <- S_design
      private$.K <- nrow(S_design)
      
      # Warn if S_design looks like it needs convolution
      if (mean(colSums(S_design) != 0) < 0.5) {
        warning("S_design appears sparse - ensure it's already convolved with HRF if using event design")
      }
      
      # Setup HRF kernel
      private$.hrf <- setup_hrf_kernel(hrf)
      if (verbose) {
        cat(sprintf("HRF kernel length: %d\n", length(private$.hrf)))
      }
      
      # Validate rank
      max_rank <- min(V, private$.K)
      if (rank > max_rank) {
        warning(sprintf("Rank %d exceeds maximum %d, setting to %d", 
                        rank, max_rank, max_rank))
        rank <- max_rank
      }
      private$.rank <- rank
      
      # Store lambda_tv
      if (lambda_tv < 0) {
        stop("lambda_tv must be non-negative")
      }
      private$.lambda_tv <- lambda_tv
      
      # Initialize cache environment
      private$.config <- list(
        cache_env = new.env(parent = emptyenv()),
        verbose = verbose
      )
      
      if (verbose) {
        cat("Learning spatial maps via GLM+SVD...\n")
      }
      
      # Learn W via GLM+SVD
      tryCatch({
        private$.learn_W_via_glm_svd(
          Y_input = private$.Y,
          S_design_input = private$.S_design,
          hrf_kernel_input = private$.hrf,
          rank_input = private$.rank
        )
      }, error = function(e) {
        if (grepl("not yet implemented", e$message)) {
          # Expected for now, set placeholder values
          if (verbose) cat("Note: GLM+SVD learning pending implementation\n")
          private$.W <- matrix(rnorm(V * private$.K), V, private$.K)
        } else {
          stop(e)
        }
      })
      
      # Pre-compute and cache WtY and WtW efficiently
      if (private$.use_low_rank || !is.null(private$.W)) {
        private$.WtY <- private$.compute_WtY()
        private$.WtW <- private$.compute_WtW()
        
        # Estimate Lipschitz constant
        private$.L_fista <- private$.estimate_lipschitz()
        
        if (verbose) {
          cat(sprintf("Lipschitz constant estimate: %.3f\n", private$.L_fista))
        }
      }
      
      # Initialize X_hat to zeros
      private$.X_hat <- matrix(0, private$.K, T)
      
      # Initialize tracking variables
      private$.objective_values <- numeric()
      private$.converged <- FALSE
      private$.iterations <- 0
      private$.fitted <- FALSE
      
      if (verbose) {
        cat("CLD initialization complete.\n")
        cat(sprintf("Data: %d voxels x %d timepoints\n", V, T))
        cat(sprintf("States: %d, Rank: %d, Lambda TV: %.4f\n", 
                    private$.K, private$.rank, private$.lambda_tv))
      }
      
      invisible(self)
    },
    
    #' @description
    #' Fit the CLD model using FISTA optimization.
    #'
    #' @param max_iter Maximum number of FISTA iterations
    #' @param tol Convergence tolerance
    #' @param verbose Logical, print iteration progress
    #'
    #' @return Self (invisibly) for method chaining
    fit = function(max_iter = 100, tol = 1e-4, verbose = FALSE) {
      # Check if initialized
      if (is.null(private$.W)) {
        stop("Model not initialized. Call $new() first.")
      }
      
      if (private$.fitted && verbose) {
        cat("Note: Model already fitted. Continuing from previous solution.\n")
      }
      
      # Store original verbose setting
      orig_verbose <- private$.config$verbose
      if (verbose) {
        private$.config$verbose <- TRUE
      }
      
      # Run FISTA optimization
      tryCatch({
        private$.fista_tv(max_iter = max_iter, tol = tol, verbose = verbose)
        
        if (verbose) {
          cat("\nFISTA optimization complete.\n")
          cat(sprintf("Converged: %s\n", if (private$.converged) "Yes" else "No"))
          cat(sprintf("Total iterations: %d\n", private$.iterations))
          
          if (length(private$.objective_values) > 0) {
            cat(sprintf("Final objective value: %.4f\n", 
                        tail(private$.objective_values, 1)))
          }
        }
      }, error = function(e) {
        # Restore verbose setting
        private$.config$verbose <- orig_verbose
        stop("Error in FISTA optimization: ", e$message)
      })
      
      # Restore verbose setting
      private$.config$verbose <- orig_verbose
      
      invisible(self)
    },
    
    #' @description
    #' Get estimated spatial maps.
    #'
    #' @param format Output format ("matrix" or "neuroim")
    #'
    #' @return Spatial maps W (V x K matrix or list of NeuroVol)
    get_spatial_maps = function(format = "matrix") {
      # Check if we have spatial maps (either low-rank or full)
      if (!private$.use_low_rank && is.null(private$.W)) {
        warning("No spatial maps available - model not initialized")
        return(NULL)
      }
      
      format <- match.arg(format, c("matrix", "neuroim"))
      
      if (format == "matrix") {
        return(private$.get_W())
      } else {  # neuroim format
        if (!is.null(private$.Y_info$space)) {
          # Convert to NeuroVol objects
          restore_spatial_structure(
            mat = private$.W,
            reference = private$.Y_info,
            output_type = "spatial"
          )
        } else {
          warning("No spatial information available - returning matrix")
          return(private$.W)
        }
      }
    },
    
    #' @description
    #' Get estimated state activations.
    #'
    #' @return State activations X_hat (K x T matrix)
    get_activations = function() {
      if (is.null(private$.X_hat)) {
        warning("No state activations available - model not fitted")
        return(NULL)
      }
      
      return(private$.X_hat)
    },
    
    #' @description
    #' Get convergence diagnostics.
    #'
    #' @return List with objective values and convergence info
    get_diagnostics = function() {
      list(
        objective_values = private$.objective_values,
        converged = private$.converged,
        iterations = private$.iterations
      )
    },
    
    #' @description
    #' Print method for CLD object.
    print = function() {
      cat("Continuous Linear Decoder (CLD)\n")
      cat("===============================\n")
      if (!is.null(private$.Y)) {
        cat("Data dimensions: ", nrow(private$.Y), " voxels x ", 
            ncol(private$.Y), " timepoints\n", sep = "")
      }
      if (!is.null(private$.K)) {
        cat("Number of states: ", private$.K, "\n", sep = "")
        cat("Rank: ", private$.rank, "\n", sep = "")
        cat("Lambda TV: ", private$.lambda_tv, "\n", sep = "")
      }
      if (!is.null(private$.W)) {
        cat("Status: Initialized", 
            if (private$.fitted) " and fitted" else "", "\n", sep = "")
      } else {
        cat("Status: Not initialized\n")
      }
      invisible(self)
    }
  ),
  
  private = list(
    # Data and parameters
    .Y = NULL,              # V x T data matrix
    .Y_info = NULL,         # Metadata from validate_fmri_input
    .S_design = NULL,       # K x T design matrix
    .hrf = NULL,            # HRF kernel vector
    .rank = NULL,           # Low-rank approximation rank
    .lambda_tv = NULL,      # TV regularization parameter
    .K = NULL,              # Number of states
    
    # Learned parameters
    .W = NULL,              # V x K spatial maps (only if not using low-rank)
    .WtY = NULL,            # K x T pre-computed W'Y for efficiency
    .WtW = NULL,            # K x K pre-computed W'W
    .X_hat = NULL,          # K x T estimated states
    .U_r = NULL,            # V x rank low-rank spatial basis
    .S_r = NULL,            # rank singular values
    .V_r = NULL,            # K x rank low-rank loadings
    .use_low_rank = FALSE   # Flag for low-rank representation
    
    # Algorithm parameters
    .L_fista = NULL,        # Lipschitz constant for FISTA
    .config = NULL,         # Configuration settings
    
    # Convergence tracking
    .objective_values = NULL,
    .converged = FALSE,
    .iterations = 0,
    .fitted = FALSE,
    
    # Private methods
    
    #' Learn spatial maps W via GLM+SVD
    .learn_W_via_glm_svd = function(Y_input, S_design_input, hrf_kernel_input, rank_input) {
      V <- nrow(Y_input)
      K <- nrow(S_design_input)
      T <- ncol(Y_input)
      
      if (private$.config$verbose) {
        cat("Step 1: Convolving design matrix with HRF...\n")
      }
      
      # Step 1: Create convolved design matrix
      X_conv <- convolve_with_hrf(S_design_input, hrf_kernel_input)
      
      # Transpose for GLM (T x K)
      X_conv_t <- t(X_conv)
      
      # Step 2: Run voxel-wise GLMs in parallel
      
      # Check OpenMP support
      openmp_info <- check_openmp_support()
      use_parallel <- openmp_info$available && V > 100  # Only parallelize for larger problems
      
      if (private$.config$verbose && use_parallel) {
        cat(sprintf("Running parallel GLM fitting with %d threads...\n", 
                    openmp_info$max_threads))
      }
      
      # Use parallel GLM fitting if available, otherwise fall back to sequential
      if (use_parallel) {
        # Use parallel C++ implementation
        B <- parallel_glm_fit_rcpp(
          Y = Y_input,
          X_conv = X_conv_t,
          n_threads = 0,  # Auto-detect
          chunk_size = max(10, V / (openmp_info$max_threads * 10))
        )
      } else {
        # Fallback to sequential implementation
        B <- matrix(0, V, K)
        
        # Check if speedglm is available
        use_speedglm <- requireNamespace("speedglm", quietly = TRUE)
        
        if (private$.config$verbose) {
          cat("Running sequential GLM fitting...\n")
          if (V > 1000) {
            pb <- txtProgressBar(min = 0, max = V, style = 3)
          }
        }
        
        for (v in seq_len(V)) {
          if (use_speedglm) {
            fit <- speedglm::speedlm.fit(
              y = Y_input[v, ],
              X = X_conv_t,
              intercept = FALSE
            )
            B[v, ] <- coef(fit)
          } else {
            fit <- lm.fit(x = X_conv_t, y = Y_input[v, ])
            B[v, ] <- coef(fit)
          }
          
          if (private$.config$verbose && V > 1000 && v %% 100 == 0) {
            setTxtProgressBar(pb, v)
          }
        }
        
        if (private$.config$verbose && V > 1000) {
          close(pb)
        }
      }
      
      # Replace any NA coefficients with 0
      B[is.na(B)] <- 0
      
      if (private$.config$verbose) {
        cat("Step 3: Performing truncated SVD...\n")
      }
      
      # Step 3: Truncated SVD with numerical stability checks
      
      # Check for numerical issues in B
      if (any(!is.finite(B))) {
        warning("Non-finite values detected in GLM coefficients, replacing with 0")
        B[!is.finite(B)] <- 0
      }
      
      # Check if B is effectively zero
      B_norm <- norm(B, "F")
      if (B_norm < .Machine$double.eps * sqrt(V * K)) {
        stop("GLM coefficient matrix is effectively zero. Check your data and design matrix.")
      }
      
      # Compute SVD with error handling
      svd_result <- tryCatch({
        svd(B, nu = rank_input, nv = rank_input)
      }, error = function(e) {
        stop("SVD failed: ", e$message, "\nThis may indicate numerical issues with the data.")
      })
      
      # Check for valid SVD output
      if (length(svd_result$d) < rank_input) {
        warning(sprintf("Matrix rank (%d) is less than requested rank (%d). Adjusting rank.",
                       length(svd_result$d), rank_input))
        rank_input <- length(svd_result$d)
      }
      
      # Extract components
      U_r <- svd_result$u[, 1:rank_input, drop = FALSE]
      S_r <- svd_result$d[1:rank_input]
      V_r <- svd_result$v[, 1:rank_input, drop = FALSE]
      
      # Check for numerical stability of singular values
      if (any(S_r < .Machine$double.eps * S_r[1])) {
        n_stable <- sum(S_r >= .Machine$double.eps * S_r[1])
        warning(sprintf("Only %d singular values are numerically stable. Consider reducing rank.",
                       n_stable))
      }
      
      # Check condition number
      condition_number <- S_r[1] / S_r[length(S_r)]
      if (condition_number > 1e10) {
        warning(sprintf("High condition number (%.2e) detected. Results may be numerically unstable.",
                       condition_number))
      }
      
      if (private$.config$verbose) {
        cat(sprintf("Singular values range: %.3f to %.3f\n", 
                    max(S_r), min(S_r)))
        var_explained <- sum(S_r^2) / sum(svd_result$d^2)
        cat(sprintf("Variance explained by rank-%d approximation: %.1f%%\n", 
                    rank_input, var_explained * 100))
      }
      
      # Step 4: Store low-rank components separately for efficiency
      # W = U_r %*% diag(S_r) %*% t(V_r)
      # We store these separately to avoid creating the full V x K matrix
      private$.U_r <- U_r
      private$.S_r <- S_r
      private$.V_r <- V_r
      
      # For compatibility, create a function to compute W on demand
      # This avoids storing the full matrix unless necessary
      private$.W <- NULL  # Clear any existing W
      
      # Set flag indicating we're using low-rank form
      private$.use_low_rank <- TRUE
      
      if (private$.config$verbose) {
        cat("Spatial map learning complete.\n")
      }
      
      invisible(TRUE)
    },
    
    #' Run FISTA optimization
    .fista_tv = function(max_iter, tol, verbose) {
      # Check prerequisites
      if (is.null(private$.WtY) || is.null(private$.W) || is.null(private$.hrf)) {
        stop("Model not properly initialized. Run initialize() first.")
      }
      
      # Call Rcpp implementation
      result <- fista_tv_rcpp(
        WtY = private$.WtY,
        W = private$.get_W(),
        hrf_kernel = private$.hrf,
        lambda_tv = private$.lambda_tv,
        L_fista = private$.L_fista,
        X_init = private$.X_hat,  # Use current X_hat as initialization
        max_iter = max_iter,
        tol = tol,
        verbose = verbose
      )
      
      # Store results
      private$.X_hat <- result$X_hat
      private$.objective_values <- c(private$.objective_values, result$objective_values)
      private$.converged <- result$converged
      private$.iterations <- private$.iterations + result$iterations
      private$.fitted <- TRUE
      
      invisible(result)
    },
    
    #' Compute FISTA gradient
    .compute_gradient = function(X_current) {
      # Call Rcpp implementation
      H_star_X <- convolve_rows_rcpp(X_current, private$.hrf)
      
      compute_gradient_fista_rcpp(
        Y_or_WtY = private$.WtY,
        W = private$.get_W(),
        H_star_X = H_star_X,
        hrf_kernel = private$.hrf,
        precomputed_WtY = TRUE
      )
    },
    
    #' Apply TV proximal operator
    .prox_tv = function(X_tilde, lambda_step) {
      # Call Rcpp implementation
      prox_tv_condat_rcpp(X_tilde, lambda_step)
    },
    
    #' Estimate Lipschitz constant
    .estimate_lipschitz = function() {
      if (private$.use_low_rank) {
        # Use stored components if available
        estimate_lipschitz_lowrank_rcpp(
          U = private$.U_r,
          S = private$.S_r,
          V = private$.V_r,
          hrf_kernel = private$.hrf
        )
      } else {
        # Use full W matrix
        estimate_lipschitz_rcpp(
          W = private$.W,
          hrf_kernel = private$.hrf
        )
      }
    },
    
    #' Get W matrix (compute from low-rank if needed)
    .get_W = function() {
      if (private$.use_low_rank) {
        # Compute W from low-rank factors
        private$.U_r %*% (private$.S_r * t(private$.V_r))
      } else {
        private$.W
      }
    },
    
    #' Compute W'Y efficiently using low-rank form
    .compute_WtY = function() {
      if (private$.use_low_rank) {
        # W'Y = V_r * S_r * U_r' * Y
        # More efficient: compute U_r' * Y first (r x T)
        UtY <- t(private$.U_r) %*% private$.Y
        # Then V_r * (S_r * UtY)
        private$.V_r %*% (private$.S_r * UtY)
      } else {
        # Standard computation
        t(private$.W) %*% private$.Y
      }
    },
    
    #' Compute W'W efficiently using low-rank form
    .compute_WtW = function() {
      if (private$.use_low_rank) {
        # W'W = V_r * S_r^2 * V_r'
        private$.V_r %*% (private$.S_r^2 * t(private$.V_r))
      } else {
        # Standard computation
        t(private$.W) %*% private$.W
      }
    }
  ),
  
  active = list(
    #' @field W Spatial maps (read-only)
    W = function() {
      private$.get_W()
    },
    
    #' @field X_hat Estimated state activations (read-only)
    X_hat = function() {
      private$.X_hat
    },
    
    #' @field hrf HRF kernel (read-only)
    hrf = function() {
      private$.hrf
    },
    
    #' @field lambda_tv TV regularization parameter
    lambda_tv = function(value) {
      if (missing(value)) {
        private$.lambda_tv
      } else {
        if (value < 0) stop("lambda_tv must be non-negative")
        private$.lambda_tv <- value
        private$.fitted <- FALSE  # Need to refit if lambda changes
      }
    },
    
    #' @field fitted Logical indicating if model has been fitted
    fitted = function() {
      private$.fitted
    }
  )
)

#' S3 method for converting CLD to list
#'
#' @param x ContinuousLinearDecoder object
#' @param ... Ignored
#' @return List representation of the CLD object
#' @export
as.list.ContinuousLinearDecoder <- function(x, ...) {
  list(
    W = x$W,
    X_hat = x$X_hat,
    hrf = x$hrf,
    lambda_tv = x$lambda_tv,
    fitted = x$fitted,
    diagnostics = x$get_diagnostics()
  )
}

#' S3 plot method for ContinuousLinearDecoder
#'
#' @param x ContinuousLinearDecoder object
#' @param type Type of plot ("convergence", "states", "maps", "all")
#' @param ... Additional arguments passed to plotting functions
#' @return NULL (plots displayed) or list of plots if type="all"
#' @export
plot.ContinuousLinearDecoder <- function(x, type = c("convergence", "states", "maps", "all"), ...) {
  type <- match.arg(type)
  
  if (type == "all") {
    # Create multi-panel plot
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    graphics::layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
    
    # 1. Convergence
    if (length(x$get_diagnostics()$objective_values) > 0) {
      plot_convergence(x$get_diagnostics()$objective_values, 
                       type = "objective", 
                       main = "FISTA Convergence", ...)
    } else {
      plot.new()
      text(0.5, 0.5, "No convergence data\n(model not fitted)", cex = 1.2)
    }
    
    # 2. Spatial maps (subset)
    W <- x$get_spatial_maps()
    if (!is.null(W)) {
      # Show up to 4 maps
      n_maps <- min(ncol(W), 4)
      plot_spatial_maps(W[, 1:n_maps, drop = FALSE], 
                        layout = c(2, 2),
                        main = "Spatial Maps", ...)
    } else {
      plot.new()
      text(0.5, 0.5, "No spatial maps available", cex = 1.2)
    }
    
    # 3. State timecourse
    X <- x$get_activations()
    if (!is.null(X) && any(X != 0)) {
      plot_state_timecourse(X, type = "line", 
                            main = "State Activations", ...)
    } else {
      plot.new()
      text(0.5, 0.5, "No state activations\n(model not fitted)", cex = 1.2)
    }
    
    invisible(NULL)
    
  } else if (type == "convergence") {
    obj_vals <- x$get_diagnostics()$objective_values
    if (length(obj_vals) == 0) {
      warning("No convergence data available - model not fitted")
      return(invisible(NULL))
    }
    
    plot_convergence(obj_vals, type = "objective", 
                     highlight_converged = if(x$get_diagnostics()$converged) 
                       x$get_diagnostics()$iterations else NULL, ...)
    
  } else if (type == "states") {
    X <- x$get_activations()
    if (is.null(X) || all(X == 0)) {
      warning("No state activations available - model not fitted")
      return(invisible(NULL))
    }
    
    plot_state_timecourse(X, ...)
    
  } else if (type == "maps") {
    W <- x$get_spatial_maps()
    if (is.null(W)) {
      warning("No spatial maps available")
      return(invisible(NULL))
    }
    
    plot_spatial_maps(W, ...)
  }
  
  invisible(NULL)
}
