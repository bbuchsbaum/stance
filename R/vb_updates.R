#' High-Level Wrapper for CBD Analysis
#'
#' @description
#' Simplified interface for running Continuous Bayesian Decoder analysis
#' on neuroimaging data with sensible defaults.
#'
#' @param Y fMRI data (matrix or neuroim2::NeuroVec)
#' @param K Number of latent states
#' @param r Spatial rank (auto-selected if NULL)
#' @param hrf_basis HRF basis specification
#' @param lambda_H_prior GMRF precision parameter
#' @param max_iter Maximum VB iterations
#' @param verbose Print progress
#' @param ... Additional arguments passed to ContinuousBayesianDecoder
#'
#' @return A fitted ContinuousBayesianDecoder object
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data <- simulate_cbd_data(V = 1000, T = 200, K = 3)
#' 
#' # Run analysis
#' result <- run_cbd_analysis(
#'   Y = data$Y,
#'   K = 3,
#'   hrf_basis = "canonical",
#'   max_iter = 50
#' )
#' 
#' # Extract results
#' spatial_maps <- result$get_spatial_maps()
#' states <- result$get_state_sequence()
#' }
#'
#' @export
run_cbd_analysis <- function(Y, K, r = NULL, hrf_basis = "canonical", 
                             lambda_H_prior = 1.0, max_iter = 100, 
                             verbose = TRUE, ...) {
  
  # Input validation
  if (missing(Y) || missing(K)) {
    stop("Both Y and K must be provided")
  }
  
  # Create and fit decoder
  decoder <- ContinuousBayesianDecoder$new(
    Y = Y,
    K = K,
    r = r,
    hrf_basis = hrf_basis,
    lambda_H_prior = lambda_H_prior,
    ...
  )
  
  # Fit model
  decoder$fit(max_iter = max_iter, verbose = verbose)
  
  return(decoder)
}

#' Variational Bayes E-Step Implementation
#'
#' @description
#' Internal function for the E-step of the VB algorithm.
#' Updates variational parameters for states and transitions.
#'
#' @param Y Data matrix (V x T)
#' @param params List of current parameter estimates
#' @param config Algorithm configuration
#'
#' @return Updated variational parameters
#'
#' @keywords internal
vb_e_step <- function(Y, params, config) {
  # Forward-backward algorithm for state inference
  result <- forward_backward_algorithm(
    Y = Y,
    U = params$U,
    V = params$V,
    H_v = params$H_v,
    hrf_basis = params$hrf_basis,
    Pi = params$Pi,
    pi0 = params$pi0,
    sigma2 = params$sigma2,
    engine = config$engine
  )
  
  return(list(
    S_gamma = result$gamma,
    S_xi = result$xi,
    log_likelihood = result$log_likelihood
  ))
}

#' Variational Bayes M-Step Implementation
#'
#' @description
#' Internal function for the M-step of the VB algorithm.
#' Updates model parameters given current variational parameters.
#'
#' @param Y Data matrix (V x T)
#' @param vb_params Current variational parameters
#' @param params Current parameter estimates
#' @param config Algorithm configuration
#'
#' @return Updated parameter estimates
#'
#' @keywords internal
vb_m_step <- function(Y, vb_params, params, config) {
  
  # Update spatial components (U, V)
  spatial_update <- update_spatial_components(
    Y = Y,
    S_gamma = vb_params$S_gamma,
    H_v = params$H_v,
    hrf_basis = params$hrf_basis,
    current_U = params$U,
    current_V = params$V,
    engine = config$engine
  )
  
  # Update HRF coefficients with GMRF prior
  H_v_update <- update_hrf_coefficients(
    Y = Y,
    S_gamma = vb_params$S_gamma,
    U = spatial_update$U,
    V = spatial_update$V,
    hrf_basis = params$hrf_basis,
    L_gmrf = params$L_gmrf,
    lambda_H_prior = params$lambda_H_prior,
    sigma2 = params$sigma2,
    engine = config$engine
  )
  
  # Update HMM parameters
  hmm_update <- update_hmm_parameters(
    S_gamma = vb_params$S_gamma,
    S_xi = vb_params$S_xi,
    prior_Pi = params$prior_Pi,
    prior_pi0 = params$prior_pi0
  )
  
  # Update noise variance
  sigma2_update <- update_noise_variance(
    Y = Y,
    S_gamma = vb_params$S_gamma,
    U = spatial_update$U,
    V = spatial_update$V,
    H_v = H_v_update,
    hrf_basis = params$hrf_basis,
    prior_sigma2 = params$sigma2_prior
  )
  
  return(list(
    U = spatial_update$U,
    V = spatial_update$V,
    H_v = H_v_update,
    Pi = hmm_update$Pi,
    pi0 = hmm_update$pi0,
    sigma2 = sigma2_update
  ))
}

#' Compute Evidence Lower BOund (ELBO)
#'
#' @description
#' Computes the ELBO for monitoring convergence of the VB algorithm.
#'
#' @param Y Data matrix (V x T)
#' @param vb_params Variational parameters
#' @param params Model parameters
#' @param config Algorithm configuration
#'
#' @return ELBO value
#'
#' @keywords internal
compute_elbo <- function(Y, vb_params, params, config) {
  
  # Data likelihood term
  likelihood_term <- compute_data_likelihood_term(
    Y = Y,
    S_gamma = vb_params$S_gamma,
    U = params$U,
    V = params$V,
    H_v = params$H_v,
    hrf_basis = params$hrf_basis,
    sigma2 = params$sigma2
  )
  
  # HMM prior terms
  hmm_prior_term <- compute_hmm_prior_term(
    S_gamma = vb_params$S_gamma,
    S_xi = vb_params$S_xi,
    Pi = params$Pi,
    pi0 = params$pi0
  )
  
  # GMRF prior term for HRF
  gmrf_prior_term <- compute_gmrf_prior_term(
    H_v = params$H_v,
    L_gmrf = params$L_gmrf,
    lambda_H_prior = params$lambda_H_prior
  )
  
  # Entropy terms
  entropy_term <- compute_entropy_term(
    S_gamma = vb_params$S_gamma,
    S_xi = vb_params$S_xi
  )
  
  elbo <- likelihood_term + hmm_prior_term + gmrf_prior_term + entropy_term
  
  return(elbo)
}

# Helper functions for parameter updates ----

#' Update Spatial Components
#' @keywords internal
update_spatial_components <- function(Y, S_gamma, H_v, hrf_basis, 
                                     current_U, current_V, engine = "cpp") {
  if (engine == "cpp") {
    # Use Rcpp implementation
    return(update_spatial_components_cpp(Y, S_gamma, H_v, hrf_basis, current_U, current_V))
  } else {
    # R implementation (placeholder)
    stop("R implementation not yet available")
  }
}

#' Update HRF Coefficients
#' @keywords internal
update_hrf_coefficients <- function(Y, S_gamma, U, V, hrf_basis, L_gmrf, 
                                   lambda_H_prior, sigma2, engine = "cpp") {
  if (engine == "cpp") {
    # Use Rcpp implementation with GMRF prior
    return(update_hrf_coefficients_gmrf_cpp(
      Y, S_gamma, U, V, hrf_basis, L_gmrf, lambda_H_prior, sigma2
    ))
  } else {
    # R implementation (placeholder)
    stop("R implementation not yet available")
  }
}

#' Update HMM Parameters
#' @keywords internal
update_hmm_parameters <- function(S_gamma, S_xi, prior_Pi = NULL, prior_pi0 = NULL) {
  K <- nrow(S_gamma)
  T_obs <- ncol(S_gamma)
  
  # Update initial state probabilities
  pi0 <- S_gamma[, 1]
  if (!is.null(prior_pi0)) {
    pi0 <- pi0 + prior_pi0 - 1  # Add Dirichlet prior
  }
  pi0 <- pi0 / sum(pi0)
  
  # Update transition matrix
  Pi <- apply(S_xi, c(1, 2), sum)
  if (!is.null(prior_Pi)) {
    Pi <- Pi + prior_Pi - 1  # Add Dirichlet prior
  }
  Pi <- Pi / rowSums(Pi)
  
  return(list(Pi = Pi, pi0 = pi0))
}

#' Update Noise Variance
#' @keywords internal
update_noise_variance <- function(Y, S_gamma, U, V, H_v, hrf_basis, prior_sigma2 = NULL) {
  # Compute expected residual sum of squares
  # This is a simplified version - full implementation would be more complex
  rss <- compute_residual_sum_squares(Y, S_gamma, U, V, H_v, hrf_basis)
  n_obs <- nrow(Y) * ncol(Y)
  
  sigma2 <- rss / n_obs
  
  if (!is.null(prior_sigma2)) {
    # Add inverse gamma prior
    sigma2 <- (rss + prior_sigma2$scale) / (n_obs + prior_sigma2$shape + 1)
  }
  
  return(sigma2)
}

# Placeholder functions for ELBO computation ----

#' @keywords internal
compute_data_likelihood_term <- function(Y, S_gamma, U, V, H_v, hrf_basis, sigma2) {
  # Compute expected log likelihood under variational distribution
  V_voxels <- nrow(Y)
  T <- ncol(Y)

  # Precompute spatial maps (V x K)
  W <- U %*% t(V)

  # Expected reconstruction for all time points
  Y_hat <- W %*% S_gamma

  # Reconstruction error
  reconstruction_error <- sum((Y - Y_hat)^2)

  # Log-likelihood (up to constant)
  log_lik <- -0.5 * V_voxels * T * log(2 * pi * sigma2)
  log_lik <- log_lik - 0.5 * reconstruction_error / sigma2

  log_lik
}

#' @keywords internal
compute_hmm_prior_term <- function(S_gamma, S_xi, Pi, pi0) {
  # Compute expected log prior for HMM using vectorized operations
  log_prior <- sum(S_gamma[, 1] * log(pi0 + 1e-10))

  if (!is.null(S_xi) && length(dim(S_xi)) == 3) {
    Pi_rep <- array(Pi + 1e-10, dim = dim(S_xi))
    log_prior <- log_prior + sum(S_xi * log(Pi_rep))
  }

  log_prior
}

#' @keywords internal
compute_gmrf_prior_term <- function(H_v, L_gmrf, lambda_H_prior) {
  # Compute GMRF prior for HRF coefficients
  # H_v is V x B matrix (voxels x basis functions)
  
  if (is.null(L_gmrf) || lambda_H_prior <= 0) {
    return(0)  # No spatial smoothing
  }
  
  # Compute quadratic form: -0.5 * lambda * sum_b H_b' * L * H_b
  log_prior <- 0
  B <- ncol(H_v)
  
  for (b in seq_len(B)) {
    h_b <- H_v[, b]
    # Quadratic form with graph Laplacian
    quadratic_form <- as.numeric(t(h_b) %*% L_gmrf %*% h_b)
    log_prior <- log_prior - 0.5 * lambda_H_prior * quadratic_form
  }
  
  return(log_prior)
}

#' @keywords internal
compute_entropy_term <- function(S_gamma, S_xi) {
  # Compute entropy of variational distribution
  entropy <- -sum(S_gamma * log(S_gamma + 1e-10))

  if (!is.null(S_xi) && length(dim(S_xi)) == 3) {
    entropy <- entropy - sum(S_xi * log(S_xi + 1e-10))
  }

  entropy
}

#' @keywords internal
compute_residual_sum_squares <- function(Y, S_gamma, U, V, H_v, hrf_basis) {

  # Compute expected residual sum of squares using matrix operations
  W <- U %*% t(V)
  Y_hat <- W %*% S_gamma
  sum((Y - Y_hat)^2)

  
  rss <- 0
  
  # Simplified computation - TODO full version would properly handle HRF convolution
  for (t in seq_len(T)) {
    y_t <- Y[, t]
    y_hat_t <- numeric(V_voxels)
    
    # Compute expected observation
    for (k in seq_len(K)) {
      gamma_kt <- S_gamma[k, t]
      if (gamma_kt > 1e-10) {
        # W_k = U %*% V[k, ]
        w_k <- U %*% V[k, ]
        y_hat_t <- y_hat_t + gamma_kt * w_k
      }
    }
    
    # Add squared residuals
    rss <- rss + sum((y_t - y_hat_t)^2)
  }
  
  return(rss)
}
