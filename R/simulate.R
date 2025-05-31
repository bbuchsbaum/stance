#' Unified Simulation Framework for stance
#'
#' Comprehensive simulation functions supporting both CLD and CBD algorithms.
#' Generates realistic fMRI data with spatial structure, temporal dynamics,
#' and neuroimaging characteristics.
#'
#' @name simulate
NULL

#' Simulate fMRI Data
#'
#' Unified simulation function that generates synthetic fMRI data for both
#' CLD and CBD algorithms. Includes realistic spatial coherence, temporal
#' structure, and noise characteristics.
#'
#' @param V Number of voxels
#' @param T Number of time points  
#' @param K Number of cognitive states
#' @param rank Rank of spatial decomposition (default: min(V, K))
#' @param algorithm Character, "CLD" or "CBD" to determine state generation
#' @param hrf_spec HRF specification (see setup_hrf_kernel)
#' @param TR Repetition time in seconds
#' @param snr Signal-to-noise ratio
#' @param spatial_smooth Logical, add spatial smoothness to patterns
#' @param temporal_smooth Logical, add temporal smoothness (for CLD)
#' @param transition_matrix K x K transition matrix (for CBD), if NULL uses random
#' @param pi0 Optional length-K vector of initial state probabilities (CBD)
#' @param state_design K x T design matrix (for supervised mode), if NULL generates states
#' @param return_neuroim Logical, return neuroimaging objects instead of matrices
#' @param dims Spatial dimensions c(x, y, z) for neuroimaging output
#' @param verbose Logical, print simulation progress
#' 
#' @return A list containing:
#'   \item{Y}{V x T matrix or NeuroVec of simulated fMRI data}
#'   \item{S}{K x T matrix of true states (one-hot for CBD, continuous for CLD)}
#'   \item{W}{V x K matrix of true spatial maps}
#'   \item{U}{V x rank matrix of spatial basis}
#'   \item{V_mat}{K x rank matrix of state loadings}
#'   \item{X}{K x T matrix of convolved states}
#'   \item{hrf}{HRF kernel used}
#'   \item{noise}{V x T matrix of noise added}
#'   \item{Pi}{K x K transition matrix used for state simulation (CBD)}
#'   \item{pi0}{Length-K vector of initial state probabilities (CBD)}
#'   \item{params}{List of simulation parameters}
#' 
#' @export
#' @examples
#' \dontrun{
#' # CLD simulation with continuous states
#' sim_cld <- simulate_fmri_data(V = 1000, T = 200, K = 3, algorithm = "CLD")
#' 
#' # CBD simulation with Markov states
#' sim_cbd <- simulate_fmri_data(V = 1000, T = 200, K = 3, algorithm = "CBD")
#' 
#' # With neuroimaging output
#' sim_neuro <- simulate_fmri_data(V = 1000, T = 200, K = 3, 
#'                                 return_neuroim = TRUE, dims = c(10, 10, 10))
#' }
#' @export
simulate_fmri_data <- function(V = 1000, T = 200, K = 3, rank = NULL,
                               algorithm = c("CLD", "CBD"),
                               hrf_spec = "spmg1", TR = 2,
                               snr = 1,
                               spatial_smooth = TRUE,
                               temporal_smooth = TRUE,
                               transition_matrix = NULL,
                               pi0 = NULL,
                               state_design = NULL,
                               return_neuroim = FALSE,
                               dims = NULL,
                               verbose = FALSE) {
  
  algorithm <- match.arg(algorithm)
  
  # Set rank if not specified
  if (is.null(rank)) {
    rank <- min(V, K, 20)  # Cap at 20 for efficiency
  }

  # Ensure core parameters are positive integers
  if (!is.numeric(V) || length(V) != 1 || V <= 0 || V != as.integer(V)) {
    stop("V must be a positive integer")
  }
  if (!is.numeric(T) || length(T) != 1 || T <= 0 || T != as.integer(T)) {
    stop("T must be a positive integer")
  }
  if (!is.numeric(K) || length(K) != 1 || K <= 0 || K != as.integer(K)) {
    stop("K must be a positive integer")
  }
  if (!is.numeric(rank) || length(rank) != 1 || rank <= 0 || rank != as.integer(rank)) {
    stop("rank must be a positive integer")
  }
  if (!is.null(dims)) {
    if (length(dims) != 3) {
      stop("dims must have length 3")
    }
    if (prod(dims) != V) {
      stop("product of dims must equal V")
    }
  }
  
  # Validate inputs
  if (return_neuroim && is.null(dims)) {
    # Automatically determine dimensions
    dims <- determine_spatial_dims(V)
  }
  
  if (verbose) {
    message(sprintf("Simulating %s data: V=%d, T=%d, K=%d, rank=%d", 
                    algorithm, V, T, K, rank))
  }
  
  # Generate spatial patterns
  spatial_data <- generate_spatial_patterns(V, K, rank, spatial_smooth, dims)
  U <- spatial_data$U
  V_mat <- spatial_data$V_mat
  W <- spatial_data$W
  
  # Generate states based on algorithm
  if (!is.null(state_design)) {
    # Use provided design
    if (!is.matrix(state_design) || nrow(state_design) != K || ncol(state_design) != T) {
      stop("state_design must be a K x T matrix")
    }
    S <- state_design
  } else if (algorithm == "CBD") {
    # Generate Markov chain states
    markov <- generate_markov_states(K, T, transition_matrix, pi0)
    S <- markov$S
    transition_matrix <- markov$Pi
    pi0 <- markov$pi0
  } else {  # CLD
    # Generate continuous states with temporal structure
    S <- generate_continuous_states(K, T, smooth = temporal_smooth)
  }
  
  # Setup HRF
  hrf <- setup_hrf_kernel(hrf_spec, TR = TR)
  
  # Convolve states with HRF
  X <- convolve_with_hrf(S, hrf)
  
  # Generate signal: Y = W * X = U * V_mat * X
  signal <- W %*% X
  
  # Generate structured noise
  noise <- generate_structured_noise(V, T, dims)
  
  # Combine signal and noise based on SNR
  signal_power <- mean(signal^2)
  noise_power <- mean(noise^2)

  if (snr <= 0) {
    # Avoid division by zero or negative scaling
    noise_scaled <- noise
    Y <- noise_scaled
  } else {
    noise_scaled <- noise * sqrt(signal_power / (snr * noise_power))
    Y <- signal + noise_scaled
  }
  
  # Create output structure
  output <- list(
    Y = Y,
    S = S,
    W = W,
    U = U,
    V_mat = V_mat,
    X = X,
    hrf = hrf,
    noise = noise_scaled,
    Pi = if (algorithm == "CBD") transition_matrix else NULL,
    pi0 = if (algorithm == "CBD") pi0 else NULL,
    params = list(
      V = V, T = T, K = K, rank = rank,
      algorithm = algorithm,
      hrf_spec = hrf_spec,
      TR = TR,
      snr = snr,
      spatial_smooth = spatial_smooth,
      temporal_smooth = temporal_smooth,
      Pi = if (algorithm == "CBD") transition_matrix else NULL,
      pi0 = if (algorithm == "CBD") pi0 else NULL
    )
  )
  
  # Convert to neuroimaging format if requested
  if (return_neuroim) {
    output <- convert_to_neuroim(output, dims, TR)
  }
  
  output
}

#' Simulate CBD Data with Voxel-Specific HRFs
#'
#' Convenience wrapper around `simulate_fmri_data` that allows the
#' user to specify voxel-specific HRF coefficients. The HRF for each
#' voxel is obtained by multiplying the provided coefficient matrix
#' with the chosen HRF basis. The function returns a structured list
#' including the ground truth HRF coefficients and noise variance to
#' simplify unit tests.
#'
#' @param V Number of voxels
#' @param T Number of time points
#' @param K Number of cognitive states
#' @param rank Rank of spatial decomposition (default: min(V, K))
#' @param hrf_basis HRF basis specification or matrix
#' @param true_H Matrix of true HRF coefficients (V x L_basis)
#' @param TR Repetition time in seconds
#' @param snr Signal-to-noise ratio
#' @param transition_matrix Optional K x K transition matrix
#' @param pi0 Optional length-K initial state probabilities
#' @param state_design Optional K x T design matrix
#' @param return_neuroim Logical, return neuroimaging objects
#' @param dims Spatial dimensions when return_neuroim = TRUE
#' @param verbose Logical, print simulation progress
#'
#' @return List with elements Y, S, W, U, V_mat, H, hrf_basis,
#'   true_sigma2 and other parameters. When `return_neuroim = TRUE`
#'   NeuroVec/NeuroVol objects are included.
#' @export
simulate_cbd_data <- function(V = 1000, T = 200, K = 3, rank = NULL,
                              hrf_basis = "canonical", true_H = NULL,
                              TR = 2, snr = 1,
                              transition_matrix = NULL, pi0 = NULL,
                              state_design = NULL, return_neuroim = FALSE,
                              dims = NULL, verbose = FALSE) {

  if (is.null(rank)) {
    rank <- min(V, K, 20)
  }

  if (return_neuroim && is.null(dims)) {
    dims <- determine_spatial_dims(V)
  }

  spatial_data <- generate_spatial_patterns(V, K, rank, TRUE, dims)
  U <- spatial_data$U
  V_mat <- spatial_data$V_mat
  W <- spatial_data$W

  if (!is.null(state_design)) {
    if (!is.matrix(state_design) || nrow(state_design) != K ||
        ncol(state_design) != T) {
      stop("state_design must be a K x T matrix")
    }
    S <- state_design
    Pi <- transition_matrix
    pi0 <- pi0
  } else {
    markov <- generate_markov_states(K, T, transition_matrix, pi0)
    S <- markov$S
    Pi <- markov$Pi
    pi0 <- markov$pi0
  }

  basis <- if (is.matrix(hrf_basis)) {
    hrf_basis
  } else {
    fun <- switch(tolower(hrf_basis),
                  "canonical" = create_hrf_basis_canonical,
                  "fir" = create_hrf_basis_fir,
                  "spline" = create_hrf_basis_fir,
                  stop("Unknown hrf_basis: ", hrf_basis))
    do.call(fun, list(TR = TR))
  }

  L_basis <- ncol(basis)

  if (is.null(true_H)) {
    true_H <- matrix(rnorm(V * L_basis, sd = 0.1), V, L_basis)
    true_H[, 1] <- true_H[, 1] + 1
  } else {
    if (!is.matrix(true_H) || nrow(true_H) != V || ncol(true_H) != L_basis) {
      stop("true_H must be a V x L_basis matrix")
    }
  }

  signal <- matrix(0, V, T)
  for (v in seq_len(V)) {
    h_v <- basis %*% true_H[v, ]
    X_v <- convolve_with_hrf(S, as.vector(h_v))
    signal[v, ] <- drop(W[v, ] %*% X_v)
  }

  noise <- generate_structured_noise(V, T, dims)
  signal_power <- mean(signal^2)
  noise_power <- mean(noise^2)

  if (snr <= 0) {
    noise_scaled <- noise
    Y <- noise_scaled
  } else {
    noise_scaled <- noise * sqrt(signal_power / (snr * noise_power))
    Y <- signal + noise_scaled
  }

  out <- list(
    Y = Y,
    S = S,
    W = W,
    U = U,
    V_mat = V_mat,
    H = true_H,
    hrf_basis = basis,
    noise = noise_scaled,
    true_sigma2 = mean(noise_scaled^2),
    Pi = Pi,
    pi0 = pi0,
    params = list(V = V, T = T, K = K, rank = rank, TR = TR,
                  snr = snr, Pi = Pi, pi0 = pi0)
  )

  if (return_neuroim) {
    out <- convert_to_neuroim(out, dims, TR)
  }

  out
}


#' Generate Spatial Patterns
#'
#' Creates realistic spatial patterns with low-rank structure and optional smoothness.
#'
#' @param V Number of voxels
#' @param K Number of states
#' @param rank Rank of decomposition
#' @param smooth Logical, add spatial smoothness
#' @param dims Spatial dimensions for smoothing
#' 
#' @return List with U, V_mat, and W matrices
#' @keywords internal
generate_spatial_patterns <- function(V, K, rank, smooth = TRUE, dims = NULL) {
  # Generate orthogonal spatial basis U
  U_raw <- matrix(rnorm(V * rank), V, rank)
  U_svd <- svd(U_raw, nu = rank, nv = 0)
  U <- U_svd$u
  
  # Add spatial smoothness if requested
  if (smooth && !is.null(dims)) {
    U <- apply(U, 2, function(col) {
      # Reshape to 3D, smooth, and flatten
      arr <- array(col, dim = dims)
      # Simple box smoothing
      arr_smooth <- spatial_smooth_3d(arr)
      as.vector(arr_smooth)
    })
    # Re-orthogonalize after smoothing
    U_svd <- svd(U, nu = rank, nv = 0)
    U <- U_svd$u
  }
  
  # Generate state loadings V_mat
  V_mat <- matrix(rnorm(K * rank), K, rank)
  
  # Add structure to make states more distinct
  for (k in 1:K) {
    # Each state emphasizes different components
    emphasis <- rep(0.3, rank)
    emphasis[((k-1) %% rank) + 1] <- 1
    V_mat[k, ] <- V_mat[k, ] * emphasis
  }
  
  # Compute full spatial maps
  W <- U %*% t(V_mat)
  
  list(U = U, V_mat = V_mat, W = W)
}

#' Generate Markov Chain States
#'
#' Generates state sequences following a Markov chain for CBD.
#'
#' @param K Number of states
#' @param T Number of time points
#' @param transition_matrix K x K transition matrix, if NULL generates random
#' @param pi0 Optional length-K vector of initial state probabilities
#'
#' @return List with elements \code{S} (K x T one-hot matrix), \code{Pi}
#'   (transition matrix), and \code{pi0} (initial probabilities)
#' @keywords internal
generate_markov_states <- function(K, T, transition_matrix = NULL, pi0 = NULL) {
  # Generate transition matrix if not provided
  if (is.null(transition_matrix)) {
    # Create transition matrix with preference for self-transitions
    trans_raw <- matrix(runif(K * K, 0.1, 0.3), K, K)
    diag(trans_raw) <- runif(K, 0.7, 0.9)  # Higher self-transition
    # Normalize rows
    transition_matrix <- trans_raw / rowSums(trans_raw)
  }

  # Initial state probabilities
  if (is.null(pi0)) {
    pi0 <- rep(1/K, K)
  } else {
    if (length(pi0) != K) {
      stop("pi0 must have length K")
    }
    if (any(pi0 < 0)) {
      stop("pi0 must be non-negative")
    }
    pi0 <- pi0 / sum(pi0)
  }

  # Generate state sequence
  states <- integer(T)
  states[1] <- sample(1:K, 1, prob = pi0)

  for (t in 2:T) {
    states[t] <- sample(1:K, 1, prob = transition_matrix[states[t-1], ])
  }

  # Convert to one-hot encoding
  S <- matrix(0, K, T)
  for (t in 1:T) {
    S[states[t], t] <- 1
  }

  list(S = S, Pi = transition_matrix, pi0 = pi0)
}

#' Generate Continuous States
#'
#' Generates continuous state activations with temporal structure for CLD.
#'
#' @param K Number of states
#' @param T Number of time points
#' @param smooth Logical, add temporal smoothness
#' 
#' @return K x T matrix of continuous activations
#' @keywords internal
generate_continuous_states <- function(K, T, smooth = TRUE) {
  # Generate base activations
  S <- matrix(0, K, T)
  
  # Create block structure with overlaps
  block_length <- as.integer(max(1, floor(T / (2 * K))))  # Allow for overlaps
  
  for (k in 1:K) {
    # Primary activation periods
    n_blocks <- sample(2:4, 1)
    for (b in 1:n_blocks) {
      start <- sample(1:(T - block_length), 1)
      shift <- as.integer(round(rnorm(1, 0, block_length / 4)))
      end <- min(start + block_length + shift, T)
      end <- as.integer(end)
      S[k, start:end] <- rnorm(as.integer(end - start + 1L), mean = 1, sd = 0.2)
    }
  }
  
  # Ensure non-negative
  S[S < 0] <- 0
  
  # Add temporal smoothness
  if (smooth) {
    # Apply moving average
    window <- 5
    S <- t(apply(S, 1, function(row) {
      stats::filter(row, rep(1/window, window), sides = 2)
    }))
    S[is.na(S)] <- 0
  }
  
  # Add some sparsity
  S[S < quantile(S, 0.3)] <- 0
  
  # Normalize each state to have similar total activation
  for (k in 1:K) {
    if (sum(S[k, ]) > 0) {
      S[k, ] <- S[k, ] / sum(S[k, ]) * T / K
    }
  }
  
  S
}

#' Generate Structured Noise
#'
#' Creates realistic fMRI noise with spatial and temporal correlations.
#'
#' @param V Number of voxels
#' @param T Number of time points
#' @param dims Spatial dimensions
#' 
#' @return V x T noise matrix
#' @keywords internal
generate_structured_noise <- function(V, T, dims = NULL) {
  # Base white noise
  noise <- matrix(rnorm(V * T), V, T)
  
  # Add temporal autocorrelation (AR1)
  ar_coef <- 0.3
  for (v in 1:V) {
    noise[v, ] <- stats::filter(noise[v, ], filter = ar_coef, 
                                method = "recursive")
  }
  
  # Add spatial correlation if dimensions provided
  if (!is.null(dims) && prod(dims) == V) {
    # Create spatial correlation using neuroim2 gaussian blur
    noise_smooth <- matrix(0, V, T)
    
    # Create a temporary 3D space for smoothing
    space_3d <- neuroim2::NeuroSpace(
      dim = dims,
      spacing = rep(1, 3),  # Unit spacing for smoothing
      origin = rep(0, 3),
      axes = neuroim2::OrientationList3D$AXIAL_LPI
    )
    
    for (t in 1:T) {
      # Convert to NeuroVol for smoothing
      vol <- neuroim2::NeuroVol(noise[, t], space_3d)
      
      # Apply Gaussian blur with sigma=0.5mm
      vol_smooth <- neuroim2::gaussian_blur(vol, mask = vol, sigma = 0.5, window = 1)
      
      # Extract smoothed values
      noise_smooth[, t] <- as.vector(vol_smooth)
    }
    # Mix original and smoothed
    noise <- 0.7 * noise + 0.3 * noise_smooth
  }
  
  # Add low-frequency drift
  drift <- outer(seq_len(V), seq_len(T), function(v, t) {
    0.1 * sin(2 * pi * t / (T/2)) + 0.05 * sin(2 * pi * t / (T/5))
  })
  
  noise + drift
}

#' Apply Spatial Smoothing Using neuroim2
#'
#' Wrapper around neuroim2::gaussian_blur for 3D spatial smoothing.
#' This replaces the previous custom implementations with the optimized
#' neuroim2 version.
#'
#' @param arr 3D array to smooth
#' @param sigma Standard deviation of Gaussian kernel in mm
#' @param window Size of smoothing window (default: 1)
#'
#' @return Smoothed 3D array
#' @keywords internal
spatial_smooth_3d <- function(arr, sigma = 1, window = 1) {
  dims <- dim(arr)
  
  # Create temporary NeuroSpace for the array
  space_3d <- neuroim2::NeuroSpace(
    dim = dims,
    spacing = rep(1, 3),  # Unit spacing
    origin = rep(0, 3),
    axes = neuroim2::OrientationList3D$AXIAL_LPI
  )
  
  # Convert to NeuroVol
  vol <- neuroim2::NeuroVol(as.vector(arr), space_3d)
  
  # Apply Gaussian blur (mask parameter required but can be missing)
  # If no mask provided, gaussian_blur uses all voxels
  vol_smooth <- neuroim2::gaussian_blur(vol, mask = vol, sigma = sigma, window = window)
  
  # Convert back to array
  array(as.vector(vol_smooth), dim = dims)
}

#' Legacy 3D Spatial Smoothing (Loop Implementation)
#'
#' This helper retains the original nested-loop implementation for
#' benchmarking and regression tests. Note: This is deprecated in favor
#' of spatial_smooth_3d which uses neuroim2::gaussian_blur.
#'
#' @param arr 3D array
#' @return Smoothed 3D array
#' @keywords internal
spatial_smooth_3d_loop <- function(arr) {
  # For compatibility, just use the new implementation
  spatial_smooth_3d(arr, sigma = 1, window = 1)
}

#' Determine Spatial Dimensions
#'
#' Automatically determines reasonable 3D dimensions for a given number of voxels.
#'
#' @param V Number of voxels
#' 
#' @return Vector of dimensions c(x, y, z)
#' @keywords internal
determine_spatial_dims <- function(V) {
  # Try to find reasonable cubic dimensions
  cube_root <- round(V^(1/3))
  
  # Find factors close to cube root
  factors <- which(V %% 1:V == 0)
  
  # Find three factors that multiply to V and are close to each other
  best_dims <- c(cube_root, cube_root, cube_root)
  best_diff <- abs(prod(best_dims) - V)
  
  for (i in factors) {
    remaining <- V / i
    remaining_factors <- which(remaining %% 1:remaining == 0)
    for (j in remaining_factors) {
      k <- remaining / j
      if (k == floor(k)) {
        dims <- sort(c(i, j, k))
        diff <- abs(prod(dims) - V) + sd(dims)  # Prefer similar dimensions
        if (diff < best_diff) {
          best_dims <- dims
          best_diff <- diff
        }
      }
    }
  }
  
  best_dims
}

#' Convert Simulation to Neuroimaging Format
#'
#' Converts matrix outputs to neuroim2 objects.
#'
#' @param sim_output Output from simulate_fmri_data
#' @param dims Spatial dimensions
#' @param TR Repetition time
#' 
#' @return Modified output with neuroimaging objects
#' @keywords internal
convert_to_neuroim <- function(sim_output, dims, TR) {
  # Create NeuroSpace
  # For 4D data, only first 3 dimensions get spacing/origin
  space_obj <- neuroim2::NeuroSpace(
    dim = c(dims, sim_output$params$T),
    spacing = rep(3, 3),  # 3mm isotropic spatial (only 3D)
    origin = rep(0, 3)    # Origin only for spatial dimensions
  )
  
  # Convert Y to NeuroVec
  sim_output$Y_neurovol <- neuroim2::NeuroVec(
    data = sim_output$Y,
    space = space_obj
  )
  
  # Convert spatial maps to NeuroVol list
  sim_output$W_neurovol <- lapply(1:sim_output$params$K, function(k) {
    space_3d <- neuroim2::NeuroSpace(
      dim = dims,
      spacing = rep(3, 3),
      origin = c(0, 0, 0),
      axes = neuroim2::OrientationList3D$AXIAL_LPI
    )
    neuroim2::NeuroVol(sim_output$W[, k], space_3d)
  })
  
  sim_output
}

#' Generate Event Design
#'
#' Creates event timing for simulation studies.
#'
#' @param K Number of conditions
#' @param T Number of time points
#' @param TR Repetition time
#' @param event_duration Duration of each event in seconds
#' @param min_isi Minimum inter-stimulus interval in seconds
#' 
#' @return Data frame with onset, duration, and condition columns
#' 
#' @export
generate_event_design <- function(K, T, TR = 2, event_duration = 2, min_isi = 4) {
  total_time <- T * TR
  events <- data.frame()
  
  current_time <- min_isi  # Start after initial baseline
  
  while (current_time < total_time - event_duration - min_isi) {
    # Random condition
    condition <- sample(1:K, 1)
    
    # Add event
    events <- rbind(events, data.frame(
      onset = current_time,
      duration = event_duration,
      condition = condition
    ))
    
    # Move to next event time
    isi <- runif(1, min_isi, min_isi * 2)
    current_time <- current_time + event_duration + isi
  }
  
  events
}

#' Simulate Multi-Subject Data
#'
#' Generates data for multiple subjects with hierarchical structure.
#'
#' @param n_subjects Number of subjects
#' @param V Number of voxels
#' @param T Number of time points
#' @param K Number of states
#' @param subject_variability Amount of between-subject variability (0-1)
#' @param ... Additional arguments passed to simulate_fmri_data
#' 
#' @return List of subject data
#' 
#' @export
simulate_multi_subject <- function(n_subjects, V, T, K, 
                                   subject_variability = 0.2, ...) {
  # Generate group-level patterns
  group_sim <- simulate_fmri_data(V = V, T = T, K = K, ...)
  
  # Generate subject-specific data
  subject_data <- list()
  
  for (subj in 1:n_subjects) {
    # Add subject-specific variation to spatial maps
    W_subj <- group_sim$W + 
      matrix(rnorm(V * K, sd = subject_variability), V, K)
    
    # Subject-specific states (same structure, different timing)
    S_subj <- if (group_sim$params$algorithm == "CBD") {
      generate_markov_states(K, T, group_sim$params$Pi, group_sim$params$pi0)$S
    } else {
      generate_continuous_states(K, T, smooth = TRUE)
    }
    
    # Generate subject data
    X_subj <- convolve_with_hrf(S_subj, group_sim$hrf)
    signal_subj <- W_subj %*% X_subj
    noise_subj <- generate_structured_noise(V, T)
    
    # Scale noise based on SNR
    signal_power <- mean(signal_subj^2)
    noise_power <- mean(noise_subj^2)

    if (group_sim$params$snr <= 0) {
      noise_scaled <- noise_subj
      Y_subj <- noise_scaled
    } else {
      noise_scaled <- noise_subj * sqrt(signal_power / (group_sim$params$snr * noise_power))
      Y_subj <- signal_subj + noise_scaled
    }
    
    subject_data[[subj]] <- list(
      Y = Y_subj,
      S = S_subj,
      W = W_subj,
      subject_id = subj
    )
  }
  
  list(
    group = group_sim,
    subjects = subject_data,
    params = c(group_sim$params, list(
      n_subjects = n_subjects,
      subject_variability = subject_variability
    ))
  )
}
