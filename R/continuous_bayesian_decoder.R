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
    #' @param use_gmrf Logical, enable spatial smoothing of HRFs
    #' @param lambda_h Precision for the GMRF prior
    #' @param mask Optional mask defining spatial neighbourhoods
    #' @param connectivity Neighbour connectivity (6, 18 or 26)
    initialize = function(Y, K, r = NULL, hrf_basis = "canonical",
                         hrf_params = list(),
                         lambda_H_prior = 1.0, sigma2_prior = 1.0,
                         engine = "cpp", use_gmrf = FALSE, lambda_h = lambda_H_prior,
                         mask = NULL, connectivity = 6) {

      # Validate basic arguments
      private$.validate_inputs(Y, K, r, hrf_basis, hrf_params)

      # Use shared utilities for data handling
      private$.setup_data_structures(Y)

      if (private$.T <= K) {
        stop("Number of timepoints T must be greater than K")
      }
      # HRF setup using shared utilities
      private$.setup_hrf_kernel(hrf_basis)
      private$.setup_hrf_basis(hrf_basis, hrf_params)

      # Store parameters
      private$.K <- K
      private$.r <- r %||% min(20, ceiling(private$.n_voxels / 10))
      private$.lambda_H_prior <- lambda_h
      private$.lambda_h <- lambda_h
      private$.sigma2_prior <- sigma2_prior
      private$.engine <- match.arg(engine, c("R", "cpp"))
      private$.use_gmrf <- isTRUE(use_gmrf)

      # Initialize model parameters
      private$.initialize_parameters()

      if (private$.use_gmrf) {
        if (is.null(mask)) {
          stop("Mask required for GMRF prior")
        }
        spatial_info <- get_spatial_neighbors(mask, connectivity)
        private$.L_gmrf <- create_gmrf_laplacian(spatial_info$neighbors,
                                                length(spatial_info$neighbors))
        Q <- private$.lambda_h * private$.L_gmrf +
             Matrix::Diagonal(nrow(private$.L_gmrf), 1e-6)
        private$.Q_chol <- Matrix::Cholesky(Q)
        private$.spatial_info <- spatial_info
      } else {
        private$.L_gmrf <- NULL
      }

      # Setup configuration
      private$.config <- list(
        max_iter = 100,
        tol = 1e-6,
        verbose = TRUE,
        save_history = TRUE,
        elbo_full = FALSE
      )

      private$.iterations <- 0
      private$.fitted <- FALSE
      
      invisible(self)
    },
    
    #' @description Fit the model using variational Bayes
    #' @param max_iter Maximum number of iterations
    #' @param tol Convergence tolerance for ELBO
    #' @param verbose Print progress information
    #' @param save_history Logical, keep full ELBO history
    fit = function(max_iter = 100, tol = 1e-4, verbose = TRUE,
                   save_history = FALSE, batch_size = NULL,
                   learning_rate = function(t) 1/(1 + 0.1 * t)) {
      private$.config$max_iter <- max_iter
      private$.config$tol <- tol
      private$.config$verbose <- verbose

      private$.elbo_history <- numeric()
      private$.converged <- FALSE
      private$.iterations <- 0

      if (verbose) {
        cat("Fitting Continuous Bayesian Decoder...\n")
        cat(sprintf("Data: %d voxels x %d timepoints\n", private$.n_voxels, private$.T))
        cat(sprintf("Model: K=%d states, r=%d rank\n", private$.K, private$.r))
        cat(sprintf("Engine: %s\n", private$.engine))
      }

      if (verbose) {
        pb <- cli::cli_progress_bar(total = max_iter)
      }

      if (is.null(batch_size)) {
        for (iter in seq_len(max_iter)) {
          log_lik <- private$.compute_log_likelihoods()
          private$.forward_backward(log_lik)

          private$.update_Pi()
          private$.update_U_V()
          private$.update_sigma2()

          elbo <- private$.compute_elbo()
          converged <- private$.check_convergence(elbo, tol)
          if (!save_history) {
            private$.elbo_history <- tail(private$.elbo_history, 2)
          }
          if (converged) {
            private$.converged <- TRUE
            private$.iterations <- iter
            break
          }

          if (verbose) cli::cli_progress_update()
          private$.iterations <- iter
        }
      } else {
        T_total <- private$.T
        if (batch_size < 100 || batch_size > T_total/2) {
          warning("Batch size should be 100-150 for optimal performance")
        }
        for (iter in seq_len(max_iter)) {
          batch_start <- sample.int(T_total - batch_size + 1L, 1L)
          batch_idx <- batch_start:(batch_start + batch_size - 1L)

          Y_batch <- private$.Y_data[, batch_idx, drop = FALSE]
          Y_proj_batch <- crossprod(private$.U, Y_batch)
          log_lik_batch <- private$.compute_log_likelihoods(Y_proj_batch)
          fb <- private$.forward_backward(log_lik_batch)

          lr <- learning_rate(iter)
          xi_stats <- apply(fb$xi, c(1, 2), sum)
          private$.Pi <- (1 - lr) * private$.Pi + lr * normalize_rows(xi_stats)

          private$.update_U_V_stochastic(Y_batch, fb$gamma, lr)
          private$.update_sigma2()

          if (iter %% 10 == 0) {
            elbo <- private$.compute_elbo()
            if (private$.check_convergence(elbo, tol)) {
              private$.converged <- TRUE
              private$.iterations <- iter
              break
            }
          }

          if (verbose) cli::cli_progress_update()
          private$.iterations <- iter
        }
      }

      if (verbose) cli::cli_progress_done()

      private$.fitted <- TRUE

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
        final_elbo = tail(private$.elbo_history, 1),
        iterations = private$.iterations
      )
    },

    #' @description
    #' Run comprehensive diagnostic checks on the fitted model.
    #' Returns a list with convergence, spatial and performance
    #' summaries.  When `output_dir` is supplied an HTML report is
    #' written to that location.
    #' @param output_dir Optional directory for HTML report
    diagnose_model_fit = function(output_dir = NULL) {
      diagnostics <- list(
        convergence = check_convergence_diagnostics(self),
        parameter_recovery = if (exists("true_params", inherits = FALSE)) {
          assess_parameter_recovery(self, true_params)
        } else {
          NULL
        },
        spatial = if (private$.use_gmrf) list(
          hrf_smoothness = assess_hrf_smoothness(self),
          effective_df = compute_effective_df_gmrf(self),
          spatial_correlation = compute_spatial_autocorrelation(self)
        ) else NULL,
        performance = list(
          iteration_time = NA_real_,
          memory_peak = NA_real_,
          meets_targets = check_performance_targets(self)
        )
      )

      if (!is.null(output_dir)) {
        generate_diagnostic_report(diagnostics, output_dir)
      }

      diagnostics
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
    .lambda_h = NULL,             # alias for precision
    .use_gmrf = FALSE,
    .Q_chol = NULL,               # Cholesky of precision matrix
    .spatial_info = NULL,         # neighbourhood metadata
    .use_batched_gmrf = FALSE,
    .sigma2_prior = NULL,         # Noise prior
    
    # Computation settings
    .engine = NULL,               # Computation engine
    .config = NULL,               # Algorithm configuration
    
    # Convergence tracking
    .converged = FALSE,
    .elbo_history = NULL,
    .iterations = 0,
    .fitted = FALSE,
    
    # Input validation
    .validate_inputs = function(Y, K, r, hrf_basis, hrf_params) {
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

      if (!is.list(hrf_params)) {
        stop("hrf_params must be a list")
      }
    },
    
    # Setup data structures
    .setup_data_structures = function(Y) {
      info <- validate_fmri_input(Y)
      extracted <- extract_data_matrix(info$data, preserve_attributes = TRUE)

      private$.Y_data <- extracted$data
      private$.neuro_metadata <- extracted$metadata
      private$.n_voxels <- nrow(private$.Y_data)
      private$.T <- ncol(private$.Y_data)
    },
    
    # Setup HRF kernel for fixed-convolution routines
    .setup_hrf_kernel = function(hrf_basis) {
      private$.hrf_kernel <- setup_hrf_kernel(hrf_basis)
    },

    # Setup HRF basis matrix
    .setup_hrf_basis = function(hrf_basis, hrf_params) {
      if (is.matrix(hrf_basis)) {
        basis <- hrf_basis
      } else {
        basis_fun <- switch(tolower(hrf_basis),
          "canonical" = create_hrf_basis_canonical,
          "fir" = create_hrf_basis_fir,
          "spline" = create_hrf_basis_fir,
          stop("Unknown hrf_basis: ", hrf_basis)
        )
        basis <- do.call(basis_fun, hrf_params)
      }
      private$.hrf_basis <- basis
      private$.hrf_kernel <- basis[, 1]
      private$.L_hrf <- nrow(basis)
      private$.L_basis <- ncol(basis)
    },
    
    # Initialize model parameters
    .initialize_parameters = function() {
      # Initialize spatial components via SVD
      Y_svd <- svd(private$.Y_data, nu = private$.r, nv = 0)
      private$.U <- Y_svd$u

      V_rand <- matrix(rnorm(private$.K * private$.r), private$.K, private$.r)
      private$.V <- qr.Q(qr(V_rand))

      private$.Y_proj <- crossprod(private$.U, private$.Y_data)
      
      # Initialize HRF coefficients with least squares
      basis <- private$.hrf_basis
      L <- min(ncol(private$.Y_data), nrow(basis))
      X <- basis[seq_len(L), , drop = FALSE]
      XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
      XtY <- crossprod(X, t(private$.Y_data[, seq_len(L), drop = FALSE]))
      private$.H_v <- t(XtX_inv %*% XtY)
      
      # Initialize HMM parameters
      private$.Pi <- diag(private$.K) * 0.8 + matrix(0.2 / private$.K, private$.K, private$.K)
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
      Y <- private$.Y_data
      U <- private$.U
      Vmat <- private$.V
      hrf <- as.vector(private$.hrf_basis %*% colMeans(private$.H_v))
      gamma <- private$.S_gamma

      lambda <- 1e-6

      # Expected design in time domain
      H_star_S <- convolve_with_hrf(gamma, hrf)

      # Low-rank representation in projected space
      X_expected <- t(Vmat) %*% H_star_S

      # Update V rows via ridge regression
      UU <- crossprod(U) + diag(lambda, ncol(U))
      UU_inv <- solve(UU)
      for (k in seq_len(private$.K)) {
        rhs <- crossprod(U, Y %*% H_star_S[k, ])
        Vmat[k, ] <- as.numeric(UU_inv %*% rhs)
      }

      # Update U using expected low-rank signal
      XX <- X_expected %*% t(X_expected) + diag(lambda, nrow(X_expected))
      U_new <- Y %*% t(X_expected) %*% solve(XX)

      # Orthonormalize columns
      U_new <- project_orthogonal(U_new)

      private$.U <- U_new
      private$.V <- Vmat
      private$.Y_proj <- crossprod(U_new, Y)
    },

    # Stochastic update of spatial components using mini-batches
    .update_U_V_stochastic = function(Y_batch, gamma_batch, lr) {
      res <- update_spatial_components_r(
        Y_batch, gamma_batch, private$.H_v, private$.hrf_basis,
        private$.U, private$.V
      )
      private$.U <- (1 - lr) * private$.U + lr * res$U
      private$.U <- project_orthogonal(private$.U)
      private$.V <- (1 - lr) * private$.V + lr * res$V
      private$.Y_proj <- crossprod(private$.U, private$.Y_data)
    },
    
    .update_hrf_coefficients = function() {
      if (isTRUE(private$.use_gmrf)) {
        private$.update_H_v_gmrf()
      } else {
        H_new <- update_hrf_coefficients(
          Y = private$.Y_data,
          S_gamma = private$.S_gamma,
          U = private$.U,
          V = private$.V,
          hrf_basis = private$.hrf_basis,
          L_gmrf = private$.L_gmrf,
          lambda_H_prior = private$.lambda_H_prior,
          sigma2 = private$.sigma2,
          engine = private$.engine
        )
        private$.H_v <- H_new
      }
    },

    .update_H_v_gmrf = function() {
      H_ls <- update_hrf_coefficients_r(
        Y = private$.Y_data,
        S_gamma = private$.S_gamma,
        U = private$.U,
        V = private$.V,
        hrf_basis = private$.hrf_basis,
        L_gmrf = NULL,
        lambda_H_prior = 0,
        sigma2 = private$.sigma2
      )
      XtX <- Matrix::Diagonal(nrow(private$.L_gmrf))
      if (exists("solve_gmrf_batched")) {
        H_sm <- solve_gmrf_batched(XtX, private$.L_gmrf, H_ls,
                                   private$.lambda_h, block_size = 64L)
        private$.H_v <- as.matrix(H_sm)
      } else {
        if (is.null(private$.Q_chol)) {
          Q <- private$.lambda_h * private$.L_gmrf +
               Matrix::Diagonal(nrow(private$.L_gmrf), 1e-6)
          private$.Q_chol <- Matrix::Cholesky(Q)
        }
        H_sm <- Matrix::solve(private$.Q_chol, H_ls)
        private$.H_v <- as.matrix(H_sm)
      }
    },
    
    .update_hmm_parameters = function() {
      gamma <- private$.S_gamma
      xi <- private$.S_xi

      pseudocount <- 1e-3
      K <- private$.K

      # Update initial distribution
      pi0 <- gamma[, 1]
      pi0 <- pi0 / sum(pi0)

      # Transition matrix
      xi_sum <- apply(xi, c(1, 2), sum)
      gamma_sum <- rowSums(gamma[, -ncol(gamma), drop = FALSE])
      Pi <- sweep(xi_sum + pseudocount, 1, gamma_sum + K * pseudocount, "/")
      Pi <- Pi / rowSums(Pi)

      private$.Pi <- Pi
      private$.pi0 <- pi0
    },

    .update_noise_variance = function() {
      W <- private$.U %*% t(private$.V)
      hrf <- as.vector(private$.hrf_basis %*% colMeans(private$.H_v))
      HX <- convolve_with_hrf(private$.S_gamma, hrf)
      Y_hat <- W %*% HX

      resid <- private$.Y_data - Y_hat
      private$.sigma2 <- mean(resid^2)
    },

    # Compute ELBO
    .compute_elbo = function() {
      log_lik <- private$.compute_log_likelihoods()

      gamma <- private$.S_gamma
      xi <- private$.S_xi
      Pi <- private$.Pi
      pi0 <- private$.pi0

      # Expected log-likelihood
      elbo_ll <- sum(gamma * log_lik)

      # Entropy terms
      entropy_gamma <- -sum(gamma * log(gamma + 1e-10))
      entropy_xi <- 0
      if (!is.null(xi) && length(dim(xi)) == 3 && dim(xi)[3] > 0) {
        entropy_xi <- -sum(xi * log(xi + 1e-10))
      }

      # HMM prior terms
      prior_term <- sum(log(pi0 + 1e-10) * gamma[, 1])
      if (!is.null(xi) && length(dim(xi)) == 3 && dim(xi)[3] > 0) {
        Pi_rep <- array(Pi + 1e-10, dim = dim(xi))
        prior_term <- prior_term + sum(xi * log(Pi_rep))
      }

      # Simple Gaussian priors on U and V
      prior_params <- -0.5 * (sum(private$.U^2) + sum(private$.V^2))

      # Spatial GMRF prior on HRF coefficients
      gmrf_term <- 0
      if (!is.null(private$.L_gmrf) && private$.lambda_H_prior > 0 &&
          !is.null(private$.H_v)) {
        for (b in seq_len(ncol(private$.H_v))) {
          h_b <- private$.H_v[, b]
          gmrf_term <- gmrf_term +
            as.numeric(t(h_b) %*% private$.L_gmrf %*% h_b)
        }

        # Add log-determinant component (approximation)
        n_vox <- nrow(private$.H_v)
        log_det <- 0.5 * (n_vox - 1) * log(private$.lambda_h)

        gmrf_prior <- log_det - 0.5 * private$.lambda_h * gmrf_term
      } else {
        gmrf_prior <- 0
      }

      hrf_prior <- -0.5 * private$.lambda_H_prior * sum(private$.H_v^2)

      elbo <- elbo_ll + entropy_gamma + entropy_xi + prior_term +
        prior_params + gmrf_prior

      if (!is.null(private$.config$elbo_full) && private$.config$elbo_full) {
        elbo <- elbo + hrf_prior
      }

      return(elbo)
    },

    # Compute log likelihoods for each state/time
    .compute_log_likelihoods = function(Y_proj_in = NULL) {
      Y_proj <- Y_proj_in %||% private$.Y_proj
      Vmat <- private$.V
      hrf <- as.vector(private$.hrf_basis %*% colMeans(private$.H_v))
      sigma2 <- private$.sigma2

      r <- private$.r
      K <- private$.K
      T_len <- ncol(Y_proj)

      if (is.null(Y_proj) || is.null(Vmat) || is.null(hrf)) {
        stop("Model parameters not initialized")
      }

      if (private$.engine == "cpp" && exists("compute_log_likelihoods_rcpp")) {
        return(compute_log_likelihoods_rcpp(Y_proj, Vmat, hrf, sigma2))
      }

      sigma2 <- max(sigma2, 1e-8)
      const_term <- -0.5 * r * log(2 * pi * sigma2)

      # Impulse responses for each time point using shared HRF utilities
      impulses <- diag(T_len)
      conv_mat <- convolve_with_hrf(impulses, hrf)
      h_at_t <- diag(conv_mat)

      log_lik <- matrix(0, K, T_len)

      for (k in seq_len(K)) { # PERF: state loop dominates likelihood time
        mu_proj_k <- outer(Vmat[k, ], h_at_t) # PERF: outer product heavy
        diff <- Y_proj - mu_proj_k
        log_lik[k, ] <- const_term - 0.5 * colSums(diff^2) / sigma2
      }

      log_lik
    },

    # Forward-backward algorithm using scaled probabilities
    .forward_backward = function(log_lik) {
      K <- nrow(log_lik)
      T_len <- ncol(log_lik)
      Pi <- private$.Pi
      pi0 <- private$.pi0

      # Validate dimensions
      if (nrow(Pi) != K || ncol(Pi) != K) {
        stop("Pi must be K x K")
      }
      if (length(pi0) != K) {
        stop("pi0 must have length K")
      }

      alpha <- matrix(0, K, T_len)
      beta <- matrix(0, K, T_len)
      c_scale <- numeric(T_len)

      if (private$.engine == "cpp" && exists("forward_pass_rcpp")) {
        fp <- forward_pass_rcpp(log_lik, Pi, pi0)
        alpha <- fp$alpha
        c_scale <- fp$c
        log_likelihood <- fp$log_likelihood
      } else {
        alpha[, 1] <- pi0 * exp(log_lik[, 1])
        c_scale[1] <- 1 / sum(alpha[, 1])
        alpha[, 1] <- alpha[, 1] * c_scale[1]

        if (T_len > 1) {
          for (t in 2:T_len) { # PERF: forward loop over time
            alpha[, t] <- (alpha[, t - 1] %*% Pi) * exp(log_lik[, t])
            c_scale[t] <- 1 / sum(alpha[, t])
            alpha[, t] <- alpha[, t] * c_scale[t]
          }
        }

        log_likelihood <- -sum(log(c_scale))
      }

      # Backward pass
      if (private$.engine == "cpp" && exists("backward_pass_rcpp")) {
        beta <- backward_pass_rcpp(log_lik, Pi, c_scale)
      } else {
        beta[, T_len] <- 1
        if (T_len > 1) {
          for (t in (T_len - 1):1) { # PERF: backward loop over time
            beta[, t] <- Pi %*% (beta[, t + 1] * exp(log_lik[, t + 1]))
            beta[, t] <- beta[, t] * c_scale[t + 1]
          }
        }
      }

      gamma <- alpha * beta

      xi <- array(0, c(K, K, max(T_len - 1, 1)))
      if (T_len > 1) {
        for (t in 1:(T_len - 1)) { # PERF: xi computation loop
          temp <- (alpha[, t] %*% t(beta[, t + 1] * exp(log_lik[, t + 1]))) * Pi
          xi[, , t] <- temp / sum(temp)
        }
      } else {
        xi <- array(numeric(0), c(K, K, 0))
      }

      private$.S_gamma <- gamma
      private$.S_xi <- xi

      return(list(gamma = gamma, xi = xi, log_likelihood = log_likelihood))
    },

    # Check convergence based on relative ELBO change
    .check_convergence = function(elbo, tol) {
      private$.elbo_history <- c(private$.elbo_history, elbo)
      n <- length(private$.elbo_history)
      if (n < 2) {
        return(FALSE)
      }
      rel_change <- abs(private$.elbo_history[n] - private$.elbo_history[n - 1]) /
        (abs(private$.elbo_history[n - 1]) + 1e-10)
      rel_change < tol
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

#' Summarize a ContinuousBayesianDecoder object
#'
#' Displays model summary statistics and average posterior state
#' probabilities.
#'
#' @param object A \code{ContinuousBayesianDecoder} instance.
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns \code{object}.
#' @export
summary.ContinuousBayesianDecoder <- function(object, ...) {
  print(object)
  cat("\nAverage state probabilities:\n")
  state_probs <- object$get_state_posteriors()
  avg_probs <- round(rowMeans(state_probs), 4)
  print(avg_probs)
  invisible(object)
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

#' Extract model coefficients from a ContinuousBayesianDecoder
#'
#' Returns the main model parameters as a named list.
#'
#' @param object A \code{ContinuousBayesianDecoder} instance.
#' @param ... Additional arguments (ignored).
#' @return Named list of parameters (U, V, Pi, pi0, sigma2, H_v).
#' @export
coef.ContinuousBayesianDecoder <- function(object, ...) {
  object$get_parameters()
}

#' Compute fitted values for a ContinuousBayesianDecoder
#'
#' Reconstructs the expected data matrix \eqn{\hat{Y} = W (h \star S_\gamma)}.
#'
#' @param object A \code{ContinuousBayesianDecoder} instance.
#' @param ... Additional arguments (ignored).
#' @return Matrix of fitted values (voxels x time).
#' @export
fitted.ContinuousBayesianDecoder <- function(object, ...) {
  W <- object$get_spatial_maps(as_neurovol = FALSE)
  states <- object$get_state_posteriors()
  hrf <- object$.__enclos_env__$private$.hrf_kernel
  H_star_S <- convolve_with_hrf(states, hrf)
  W %*% H_star_S
}
