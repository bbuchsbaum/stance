# Additional Edge Case Tests

test_that("CLD handles boundary conditions", {
  # Single voxel
  Y_single <- matrix(rnorm(1 * 100), 1, 100)
  S_single <- matrix(rnorm(3 * 100), 3, 100)
  
  cld_single <- ContinuousLinearDecoder$new(
    Y = Y_single,
    S_design = S_single,
    verbose = FALSE
  )
  expect_s3_class(cld_single, "ContinuousLinearDecoder")
  expect_equal(nrow(cld_single$W), 1)
  
  # Single state
  Y_one_state <- matrix(rnorm(100 * 50), 100, 50)
  S_one_state <- matrix(1, 1, 50)
  
  cld_one <- ContinuousLinearDecoder$new(
    Y = Y_one_state,
    S_design = S_one_state,
    verbose = FALSE
  )
  expect_s3_class(cld_one, "ContinuousLinearDecoder")
  expect_equal(ncol(cld_one$W), 1)
  
  # Very short time series
  Y_short <- matrix(rnorm(100 * 10), 100, 10)
  S_short <- matrix(rnorm(3 * 10), 3, 10)
  
  cld_short <- ContinuousLinearDecoder$new(
    Y = Y_short,
    S_design = S_short,
    verbose = FALSE
  )
  expect_s3_class(cld_short, "ContinuousLinearDecoder")
})

test_that("CLD handles special matrix conditions", {
  # Orthogonal design matrix
  S_orth <- qr.Q(qr(matrix(rnorm(100 * 3), 100, 3)))
  S_orth <- t(S_orth)  # K x T
  Y <- matrix(rnorm(200 * 100), 200, 100)
  
  cld_orth <- ContinuousLinearDecoder$new(
    Y = Y,
    S_design = S_orth,
    verbose = FALSE
  )
  expect_s3_class(cld_orth, "ContinuousLinearDecoder")
  
  # Sparse design matrix
  S_sparse <- matrix(0, 3, 100)
  S_sparse[1, 1:30] <- 1
  S_sparse[2, 31:60] <- 1
  S_sparse[3, 61:90] <- 1
  
  cld_sparse <- ContinuousLinearDecoder$new(
    Y = Y,
    S_design = S_sparse,
    verbose = FALSE
  )
  expect_s3_class(cld_sparse, "ContinuousLinearDecoder")
  
  # Design with negative values
  S_negative <- matrix(rnorm(3 * 100), 3, 100)
  cld_negative <- ContinuousLinearDecoder$new(
    Y = Y,
    S_design = S_negative,
    verbose = FALSE
  )
  expect_s3_class(cld_negative, "ContinuousLinearDecoder")
})

test_that("HRF utilities handle edge cases", {
  # Zero HRF
  hrf_zero <- rep(0, 16)
  expect_error(setup_hrf_kernel(hrf_zero), "is all zeros")
  
  # Single point HRF
  hrf_single <- 1
  result_single <- setup_hrf_kernel(hrf_single, normalize = FALSE)
  expect_equal(length(result_single), 1)
  
  # Very long HRF
  hrf_long <- rnorm(100)
  result_long <- setup_hrf_kernel(hrf_long)
  expect_equal(length(result_long), 100)
  expect_equal(sum(result_long^2), 1, tolerance = 1e-10)
})

test_that("Visualization functions handle edge cases", {
  # Empty data
  expect_warning(plot_convergence(numeric(0)), "No convergence values")
  
  # Single state
  W_single <- matrix(rnorm(100), 100, 1)
  expect_silent({
    pdf(NULL)
    plot_spatial_maps(W_single)
    dev.off()
  })
  
  # Many states
  W_many <- matrix(rnorm(100 * 20), 100, 20)
  expect_silent({
    pdf(NULL)
    plot_spatial_maps(W_many)  # Should handle layout automatically
    dev.off()
  })
  
  # Single time point
  X_single_t <- matrix(1, 3, 1)
  expect_silent({
    pdf(NULL)
    plot_state_timecourse(X_single_t)
    dev.off()
  })
})

test_that("Simulation handles extreme parameters", {
  # Zero voxels (should error)
  expect_error(simulate_fmri_data(V = 0, T = 100, K = 3))
  
  # Zero time points (should error)
  expect_error(simulate_fmri_data(V = 100, T = 0, K = 3))
  
  # Zero states (should error)
  expect_error(simulate_fmri_data(V = 100, T = 100, K = 0))
  
  # Very high rank
  sim_high_rank <- simulate_fmri_data(V = 100, T = 50, K = 3, rank = 50)
  expect_equal(ncol(sim_high_rank$U), min(50, 100, 3))
  
  # Zero SNR (pure noise)
  sim_zero_snr <- simulate_fmri_data(V = 100, T = 50, K = 3, snr = 0)
  # Signal should be dominated by noise and finite
  signal_var <- var(as.vector(sim_zero_snr$W %*% sim_zero_snr$X))
  noise_var <- var(as.vector(sim_zero_snr$noise))
  expect_true(noise_var > signal_var * 10)  # Noise much larger
  expect_true(all(is.finite(sim_zero_snr$Y)))
})

test_that("Memory efficiency for large problems", {
  # Test that pre-computed values are used efficiently
  sim <- simulate_fmri_data(V = 1000, T = 100, K = 3)
  
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    verbose = FALSE
  )
  
  # Check that pre-computed values exist
  expect_true(!is.null(cld$.__enclos_env__$private$.WtY))
  expect_true(!is.null(cld$.__enclos_env__$private$.WtW))
  
  # Dimensions should be reduced
  expect_equal(dim(cld$.__enclos_env__$private$.WtY), c(3, 100))  # K x T
  expect_equal(dim(cld$.__enclos_env__$private$.WtW), c(3, 3))     # K x K
})

test_that("CLD handles various HRF specifications", {
  sim <- simulate_fmri_data(V = 100, T = 50, K = 2)
  
  # String specification
  cld1 <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    hrf = "spmg1",
    verbose = FALSE
  )
  expect_s3_class(cld1, "ContinuousLinearDecoder")
  
  # Numeric vector
  hrf_custom <- c(0, 0.2, 0.5, 0.8, 0.6, 0.3, 0.1)
  cld2 <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    hrf = hrf_custom,
    verbose = FALSE
  )
  expect_s3_class(cld2, "ContinuousLinearDecoder")
  
  # Different basis functions
  for (basis in c("spmg2", "gamma", "gaussian")) {
    cld_basis <- ContinuousLinearDecoder$new(
      Y = sim$Y,
      S_design = sim$S,
      hrf = basis,
      verbose = FALSE
    )
    expect_s3_class(cld_basis, "ContinuousLinearDecoder")
  }
})

test_that("Error messages are informative", {
  Y <- matrix(rnorm(100 * 50), 100, 50)
  S <- matrix(rnorm(3 * 50), 3, 50)
  
  # Test various error conditions with clear messages
  expect_error(
    ContinuousLinearDecoder$new(Y = "not a matrix", S_design = S),
    "must be a matrix"
  )
  
  expect_error(
    ContinuousLinearDecoder$new(Y = Y, S_design = matrix(1, 3, 30)),
    "30 columns but Y has 50 timepoints"
  )
  
  expect_error(
    ContinuousLinearDecoder$new(Y = Y, S_design = S, lambda_tv = -1),
    "lambda_tv must be non-negative"
  )
  
  expect_error(
    setup_hrf_kernel("invalid_hrf_type"),
    "Unknown HRF specification"
  )
})