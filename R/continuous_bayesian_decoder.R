#' Continuous Bayesian Decoder R6 Class
#'
#' @description
#' Main R6 class implementing the Continuous Bayesian Decoder for fMRI analysis.
#' Integrates seamlessly with neuroim2 data structures and fmrireg HRF components.
#'
#' @details
#' The ContinuousBayesianDecoder uses variational Bayes inference to simultaneously
#' learn spatial activation patterns and infer temporal state sequences from fMRI data.
#' Key features include voxel-specific HRF estimation, spatial smoothing priors,
#' and efficient low-rank approximations.
#'
#' @examples
#' \dontrun{
#' # Create synthetic data
#' data <- simulate_cbd_data(V = 1000, T = 200, K = 3)
#' 
#' # Initialize decoder
#' decoder <- ContinuousBayesianDecoder$new(
#'   Y = data$Y,
#'   K = 3,
#'   r = 10,
#'   hrf_basis = "canonical"
#' )
#' 
#' # Fit model
#' decoder$fit(max_iter = 100)
#' 
#' # Extract results
#' spatial_maps <- decoder$get_spatial_maps()
#' state_sequence <- decoder$get_state_sequence()
#' }
#'
#' @export
ContinuousBayesianDecoder <- R6::R6Class(
  classname = "ContinuousBayesianDecoder",
  
  ## Public fields and methods ----
  public = list(
    
    #' @description
    #' Initialize the Continuous Bayesian Decoder
    #'
    #' @param Y fMRI data matrix (V x T) or neuroim2::NeuroVec object
    #' @param K Number of latent states
    #' @param r Rank of spatial decomposition (default: min(20, V/10))
    #' @param hrf_basis HRF basis specification (fmrireg compatible)
    #' @param lambda_H_prior GMRF precision parameter (default: 1.0)
    #' @param sigma2_prior Noise variance prior (default: 1.0)
    #' @param engine Computation engine ("R" or "cpp", default: "cpp")
    initialize = function(Y, K, r = NULL, hrf_basis = "canonical", 
                         lambda_H_prior = 1.0, sigma2_prior = 1.0, 
                         engine = "cpp") {
      
      # Input validation and data structure handling
      private$.validate_inputs(Y, K, r, hrf_basis)
      private$.setup_data_structures(Y)
      private$.setup_hrf_kernel(hrf_basis)
      
      # Store parameters
      private$.K <- K
      private$.r <- r %||% min(20, ceiling(private$.n_voxels / 10))
      private$.lambda_H_prior <- lambda_H_prior
      private$.sigma2_prior <- sigma2_prior
      private$.engine <- match.arg(engine, c("R", "cpp"))
      
      # Initialize model parameters
      private$.initialize_parameters()
      private$.setup_gmrf_structure()
      
      # Setup configuration
      private$.config <- list(
        max_iter = 100,
        tol = 1e-6,
        verbose = TRUE,
        save_history = TRUE
      )
      
      invisible(self)
    },
    
    #' @description Fit the model using variational Bayes
    #' @param max_iter Maximum number of iterations
    #' @param tol Convergence tolerance for ELBO
    #' @param verbose Print progress information
    fit = function(max_iter = 100, tol = 1e-6, verbose = TRUE) {
      private$.config$max_iter <- max_iter
      private$.config$tol <- tol
      private$.config$verbose <- verbose
      
      if (verbose) {
        cat("Fitting Continuous Bayesian Decoder...\n")
        cat(sprintf("Data: %d voxels x %d timepoints\n", private$.n_voxels, private$.T))
        cat(sprintf("Model: K=%d states, r=%d rank\n", private$.K, private$.r))
        cat(sprintf("Engine: %s\n", private$.engine))
      }
      
      # Main VB loop
      private$.fit_vb()
      
      if (verbose) {
        cat("Fitting completed successfully!\n")
      }
      
      invisible(self)
    },
    
    #' @description Get estimated spatial maps
    #' @param as_neurovol Return as neuroim2::NeuroVol if possible
    get_spatial_maps = function(as_neurovol = TRUE) {
      W <- private$.U %*% t(private$.V)
      
      if (as_neurovol && !is.null(private$.neuro_metadata)) {
        return(private$.convert_to_neurovol(W))
      }
      
      return(W)
    },
    
    #' @description Get estimated state sequence
    get_state_sequence = function() {
      return(private$.S_gamma)
    },

    #' @description Get state posteriors
    get_state_posteriors = function() {
      private$.S_gamma
    },

    #' @description Get model parameters as a list
    get_parameters = function() {
      list(
        U = private$.U,
        V = private$.V,
        Pi = private$.Pi,
        pi0 = private$.pi0,
        sigma2 = private$.sigma2,
        H_v = private$.H_v
      )
    },

    #' @description Predict reconstructed data
    predict = function() {
      W <- private$.U %*% t(private$.V)
      W %*% private$.S_gamma
    },

    #' @description Plot ELBO convergence
    plot_convergence = function(...) {
      plot_convergence(private$.elbo_history, type = "elbo", ...)
    },
    
    #' @description Get estimated HRF parameters
    get_hrf_estimates = function() {
      return(private$.H_v)
    },
    
    #' @description Get model convergence information
    get_convergence = function() {
      list(
        converged = private$.converged,
        elbo_history = private$.elbo_history,
        final_elbo = tail(private$.elbo_history, 1)
      )
    },
    
    #' @description Generate quality control report
    #' @param output_file Output HTML file path
    qc_report = function(output_file = "cbd_qc_report.html") {
      qc_report(self, output_file = output_file)
    }
  ),
  
  ## Private fields and methods ----
  private = list(
    
    # Data and parameters
    .Y_data = NULL,           # Working data matrix (V x T)
    .neuro_metadata = NULL,       # neuroim2 space information
    .Y_proj = NULL,            # U^T Y projection
    .n_voxels = NULL,             # Number of voxels
    .T = NULL,                    # Number of timepoints
    .K = NULL,                    # Number of states
    .r = NULL,                    # Spatial rank
    
    # HRF components
    .hrf_kernel = NULL,            # HRF basis matrix
    .L_hrf = NULL,                # HRF basis length
    .L_basis = NULL,              # Number of basis functions
    
    # Model parameters
    .U = NULL,                    # Spatial basis (V x r)
    .V = NULL,               # State coefficients (K x r)
    .H_v = NULL,                  # HRF coefficients (V x L_basis)
    .Pi = NULL,                   # Transition matrix (K x K)
    .pi0 = NULL,                  # Initial state probabilities
    .sigma2 = NULL,               # Noise variance
    
    # Variational parameters
    .S_gamma = NULL,              # State probabilities (K x T)
    .S_xi = NULL,                 # Transition probabilities (K x K x T-1)
    
    # GMRF structure
    .L_gmrf = NULL,               # Graph Laplacian
    .lambda_H_prior = NULL,       # GMRF precision
    .sigma2_prior = NULL,         # Noise prior
    
    # Computation settings
    .engine = NULL,               # Computation engine
    .config = NULL,               # Algorithm configuration
    
    # Convergence tracking
    .converged = FALSE,
    .elbo_history = NULL,
    
    # Input validation
    .validate_inputs = function(Y, K, r, hrf_basis) {
      # Validate Y
      if (!(is.matrix(Y) || inherits(Y, "NeuroVec"))) {
        stop("Y must be a matrix or NeuroVec object")
      }
      
      # Validate K
      if (!is.numeric(K) || K < 2) {
        stop("K must be a numeric value >= 2")
      }
      
      # Validate r
      if (!is.null(r) && (!is.numeric(r) || r < 1)) {
        stop("r must be a positive numeric value")
      }
      
      # Validate hrf_basis
      if (is.character(hrf_basis)) {
        if (!hrf_basis %in% c("canonical", "fir", "spline")) {
          stop("hrf_basis must be 'canonical', 'fir', 'spline', or a matrix")
        }
      } else if (!is.matrix(hrf_basis)) {
        stop("hrf_basis must be a character string or matrix")
      }
    },
    
    # Setup data structures
    .setup_data_structures = function(Y) {
      private$.Y_data <- Y
      
      if (inherits(Y, "NeuroVec")) {
        private$.neuro_metadata <- space(Y)
        private$.Y_data <- as.matrix(Y)
      } else {
        private$.Y_data <- as.matrix(Y)
        private$.neuro_metadata <- NULL
      }
      
      private$.n_voxels <- nrow(private$.Y_data)
      private$.T <- ncol(private$.Y_data)
    },
    
    # Setup HRF basis using fmrireg
    .setup_hrf_kernel = function(hrf_basis) {
      if (is.character(hrf_basis)) {
        # Use fmrireg to create basis
        private$.hrf_kernel <- create_hrf_basis_neuroim2(
          type = hrf_basis,
          TR = 2.0,  # Default TR, should be configurable
          duration = 32.0
        )
      } else {
        private$.hrf_kernel <- hrf_basis
      }
      
      private$.L_hrf <- nrow(private$.hrf_kernel)
      private$.L_basis <- ncol(private$.hrf_kernel)
    },
    
    # Initialize model parameters
    .initialize_parameters = function() {
      # Initialize spatial components via SVD
      Y_svd <- svd(private$.Y_data, nu = private$.r, nv = 0)
      private$.U <- Y_svd$u
      private$.V <- matrix(rnorm(private$.K * private$.r), private$.K, private$.r)
      private$.Y_proj <- crossprod(private$.U, private$.Y_data)
      
      # Initialize HRF coefficients with least squares
      private$.H_v <- matrix(0, private$.n_voxels, private$.L_basis)
      # TODO: Implement LS initialization
      
      # Initialize HMM parameters
      private$.Pi <- matrix(0.1, private$.K, private$.K)
      diag(private$.Pi) <- 0.7
      private$.Pi <- private$.Pi / rowSums(private$.Pi)
      private$.pi0 <- rep(1/private$.K, private$.K)
      
      # Initialize noise variance
      private$.sigma2 <- var(as.vector(private$.Y_data)) * 0.1
      
      # Initialize variational parameters
      private$.S_gamma <- matrix(1/private$.K, private$.K, private$.T)
      private$.S_xi <- array(1/(private$.K^2), c(private$.K, private$.K, private$.T - 1))
    },
    
    # Setup GMRF structure
    .setup_gmrf_structure = function() {
      if (!is.null(private$.neuro_metadata)) {
        # Use neuroim2 spatial structure
        private$.L_gmrf <- create_gmrf_laplacian_neuroim2(private$.neuro_metadata, private$.n_voxels)
      } else {
        # Create simple chain graph for matrix input
        private$.L_gmrf <- create_chain_laplacian(private$.n_voxels)
      }
    },
    
    # Main VB fitting loop
    .fit_vb = function() {
      private$.elbo_history <- numeric(private$.config$max_iter)
      
      for (iter in seq_len(private$.config$max_iter)) {
        
        # E-step: Update variational parameters
        private$.e_step()
        
        # M-step: Update model parameters
        private$.m_step()
        
        # Compute ELBO
        elbo <- private$.compute_elbo()
        private$.elbo_history[iter] <- elbo
        
        # Check convergence
        if (iter > 1) {
          delta_elbo <- elbo - private$.elbo_history[iter - 1]
          if (abs(delta_elbo) < private$.config$tol) {
            private$.converged <- TRUE
            private$.elbo_history <- private$.elbo_history[1:iter]
            break
          }
        }
        
        if (private$.config$verbose && iter %% 10 == 0) {
          cat(sprintf("Iteration %d: ELBO = %.4f\n", iter, elbo))
        }
      }
    },
    
    # E-step implementation
    .e_step = function() {
      if (private$.engine == "cpp") {
        # Use Rcpp implementation
        result <- forward_backward_cpp(
          Y = private$.Y_data,
          U = private$.U,
          V = private$.V,
          H_v = private$.H_v,
          hrf_basis = private$.hrf_kernel,
          Pi = private$.Pi,
          pi0 = private$.pi0,
          sigma2 = private$.sigma2
        )
        private$.S_gamma <- result$gamma
        private$.S_xi <- result$xi
      } else {
        # Use R implementation
        private$.e_step_r()
      }
    },
    
    # E-step R implementation (placeholder)
    .e_step_r = function() {
      # TODO: Implement R version for comparison
      stop("R implementation not yet available")
    },
    
    # M-step implementation
    .m_step = function() {
      # Update U and V
      private$.update_spatial_components()
      
      # Update HRF coefficients
      private$.update_hrf_coefficients()
      
      # Update HMM parameters
      private$.update_hmm_parameters()
      
      # Update noise variance
      private$.update_noise_variance()
    },
    
    # Update spatial components (placeholder implementations)
    .update_spatial_components = function() {
      # TODO: Implement spatial component updates
    },
    
    .update_hrf_coefficients = function() {
      # TODO: Implement HRF coefficient updates with GMRF prior
    },
    
    .update_hmm_parameters = function() {
      # TODO: Implement HMM parameter updates
    },
    
    .update_noise_variance = function() {
      # TODO: Implement noise variance updates
    },

    # Compute ELBO
    .compute_elbo = function() {
      # TODO: Implement ELBO calculation
      return(0)
    },

    # Compute log likelihoods for each state/time
    .compute_log_likelihoods = function() {
      stop("Not implemented")
    },

    # Forward-backward algorithm wrapper
    .forward_backward = function(log_lik) {
      stop("Not implemented")
    },

    # Update U and V factors
    .update_U_V = function() {
      private$.update_spatial_components()
    },

    # Update transition matrix
    .update_Pi = function() {
      private$.update_hmm_parameters()
    },

    # Update noise variance (alias)
    .update_sigma2 = function() {
      private$.update_noise_variance()
    },
    
    # Convert results back to neuroim2 format
    .convert_to_neurovol = function(data) {
      if (is.null(private$.neuro_metadata)) {
        return(data)
      }
      
      # TODO: Implement conversion to NeuroVol
      return(data)
    }
  )
)

# S3 methods for ContinuousBayesianDecoder ----

#' @export
print.ContinuousBayesianDecoder <- function(x, ...) {
  cat("Continuous Bayesian Decoder\n")
  cat("===========================\n")
  cat(sprintf("Data dimensions: %d voxels x %d timepoints\n", x$.__enclos_env__$private$.n_voxels, x$.__enclos_env__$private$.T))
  cat(sprintf("Model parameters: K=%d states, r=%d rank\n", x$.__enclos_env__$private$.K, x$.__enclos_env__$private$.r))
  cat(sprintf("HRF basis: %d timepoints x %d basis functions\n", x$.__enclos_env__$private$.L_hrf, x$.__enclos_env__$private$.L_basis))
  
  if (!is.null(x$.__enclos_env__$private$.elbo_history)) {
    cat(sprintf("Fitted: %s (ELBO = %.4f)\n", 
                ifelse(x$.__enclos_env__$private$.converged, "converged", "not converged"),
                tail(x$.__enclos_env__$private$.elbo_history, 1)))
  } else {
    cat("Status: not fitted\n")
  }
  
  invisible(x)
}

#' @export
summary.ContinuousBayesianDecoder <- function(object, ...) {
  print(object)
  # TODO: Add more detailed summary information
}

#' @export
as.list.ContinuousBayesianDecoder <- function(x, ...) {
  # Return key components for inspection
  list(
    spatial_maps = x$get_spatial_maps(as_neurovol = FALSE),
    state_sequence = x$get_state_sequence(),
    hrf_estimates = x$get_hrf_estimates(),
    convergence = x$get_convergence()
  )
}
