# Robustness Testing and Edge Cases (S1a-T18)

test_that("Condat TV prox handles edge cases correctly", {
  skip_if_not(exists("prox_tv_condat_rcpp"), "Rcpp functions not compiled")
  
  # Test 1: Constant signal should be unchanged
  x_const <- rep(5, 100)
  result_const <- stance:::prox_tv_condat_rcpp(matrix(x_const, nrow = 1), lambda_tv = 0.5)
  expect_equal(as.vector(result_const), x_const, tolerance = 1e-10)
  
  # Test 2: Linear ramp should be unchanged for small lambda
  x_ramp <- seq(0, 10, length.out = 100)
  result_ramp <- stance:::prox_tv_condat_rcpp(matrix(x_ramp, nrow = 1), lambda_tv = 0.001)
  expect_equal(as.vector(result_ramp), x_ramp, tolerance = 0.01)
  
  # Test 3: Step function should be smoothed
  x_step <- c(rep(0, 50), rep(1, 50))
  result_step <- stance:::prox_tv_condat_rcpp(matrix(x_step, nrow = 1), lambda_tv = 0.2)
  # Should have fewer jumps than original
  tv_original <- sum(abs(diff(x_step)))
  tv_result <- sum(abs(diff(as.vector(result_step))))
  expect_true(tv_result < tv_original)
  
  # Test 4: Impulse noise on smooth signal should be removed
  x_smooth <- sin(seq(0, 4*pi, length.out = 100))
  x_noisy <- x_smooth
  x_noisy[c(25, 50, 75)] <- x_noisy[c(25, 50, 75)] + 2  # Add spikes
  result_denoised <- stance:::prox_tv_condat_rcpp(matrix(x_noisy, nrow = 1), lambda_tv = 0.5)
  # Should be closer to smooth signal
  error_noisy <- mean((x_noisy - x_smooth)^2)
  error_denoised <- mean((as.vector(result_denoised) - x_smooth)^2)
  expect_true(error_denoised < error_noisy)
})

test_that("CLD handles numerical stability issues", {
  # Test with very high SNR data
  set.seed(123)
  sim_high_snr <- simulate_fmri_data(V = 100, T = 50, K = 2, snr = 1000)
  cld_high <- ContinuousLinearDecoder$new(
    Y = sim_high_snr$Y,
    S_design = sim_high_snr$S,
    verbose = FALSE
  )
  expect_s3_class(cld_high, "ContinuousLinearDecoder")
  
  # Test with very low SNR data
  sim_low_snr <- simulate_fmri_data(V = 100, T = 50, K = 2, snr = 0.01)
  cld_low <- ContinuousLinearDecoder$new(
    Y = sim_low_snr$Y,
    S_design = sim_low_snr$S,
    verbose = FALSE
  )
  expect_s3_class(cld_low, "ContinuousLinearDecoder")
  
  # Test with poorly conditioned W
  V <- 100
  K <- 5
  # Create nearly collinear columns
  W_poor <- matrix(rnorm(V * K), V, K)
  W_poor[, 2] <- W_poor[, 1] + rnorm(V, sd = 0.001)  # Almost identical to column 1
  
  S_design <- matrix(rnorm(K * 50), K, 50)
  Y <- W_poor %*% S_design + matrix(rnorm(V * 50, sd = 0.1), V, 50)
  
  cld_poor <- ContinuousLinearDecoder$new(
    Y = Y,
    S_design = S_design,
    verbose = FALSE
  )
  expect_s3_class(cld_poor, "ContinuousLinearDecoder")
  
  # Test with long time series
  sim_long <- simulate_fmri_data(V = 50, T = 5000, K = 2)
  cld_long <- ContinuousLinearDecoder$new(
    Y = sim_long$Y,
    S_design = sim_long$S,
    verbose = FALSE
  )
  expect_s3_class(cld_long, "ContinuousLinearDecoder")
})

test_that("Data structure conversion preserves information", {
  # Test conversion between matrix and neuroim2 objects
  V <- 1000
  T <- 100
  dims <- c(10, 10, 10)
  
  # Create test data
  Y_mat <- matrix(rnorm(V * T), V, T)
  
  # Test with spatial metadata
  metadata <- list(
    space = list(dim = dims, spacing = c(3, 3, 3, 2)),
    mask = sample(c(TRUE, FALSE), V, replace = TRUE, prob = c(0.8, 0.2))
  )
  
  # Extract and restore
  extracted <- extract_data_matrix(Y_mat)
  expect_equal(extracted$data, Y_mat)
  
  # Test mask handling
  Y_masked <- Y_mat
  Y_masked[!metadata$mask, ] <- NA
  validated <- validate_fmri_input(Y_masked, check_finite = FALSE)
  expect_equal(dim(validated$data), dim(Y_mat))
})

test_that("CLD handles missing data appropriately", {
  # Create data with NaN values
  Y <- matrix(rnorm(100 * 50), 100, 50)
  Y[sample(length(Y), 50)] <- NaN
  S <- matrix(rnorm(3 * 50), 3, 50)
  
  # Should error with check_finite = TRUE (default)
  expect_error(
    ContinuousLinearDecoder$new(Y = Y, S_design = S),
    "non-finite values"
  )
  
  # Test with all zeros in a voxel (common in masked data)
  Y_zeros <- matrix(rnorm(100 * 50), 100, 50)
  Y_zeros[1:10, ] <- 0  # First 10 voxels are all zeros
  
  cld_zeros <- ContinuousLinearDecoder$new(
    Y = Y_zeros,
    S_design = S,
    verbose = FALSE
  )
  expect_s3_class(cld_zeros, "ContinuousLinearDecoder")
})

test_that("CLD input validation catches errors", {
  # Test mismatched dimensions
  Y <- matrix(rnorm(100 * 50), 100, 50)
  
  # Wrong K dimension
  S_wrong_K <- matrix(rnorm(3 * 40), 3, 40)
  expect_error(
    ContinuousLinearDecoder$new(Y = Y, S_design = S_wrong_K),
    "S_design has 40 columns but Y has 50 timepoints"
  )
  
  # Non-matrix S_design
  expect_error(
    ContinuousLinearDecoder$new(Y = Y, S_design = list(a = 1)),
    "S_design must be a matrix"
  )
  
  # Invalid rank
  S <- matrix(rnorm(3 * 50), 3, 50)
  expect_warning(
    ContinuousLinearDecoder$new(Y = Y, S_design = S, rank = 200),
    "Rank 200 exceeds maximum"
  )
  
  # Negative lambda
  expect_error(
    ContinuousLinearDecoder$new(Y = Y, S_design = S, lambda_tv = -0.1),
    "lambda_tv must be non-negative"
  )
})

test_that("HRF convolution handles edge cases", {
  # Empty signal
  X_empty <- matrix(nrow = 0, ncol = 0)
  hrf <- setup_hrf_kernel("spmg1")
  expect_error(convolve_with_hrf(X_empty, hrf))
  
  # Single time point
  X_single <- matrix(1, nrow = 2, ncol = 1)
  result_single <- convolve_with_hrf(X_single, hrf)
  expect_equal(dim(result_single), c(2, 1))
  
  # Very long HRF
  hrf_long <- rep(0.1, 100)
  X <- matrix(rnorm(3 * 50), 3, 50)
  result_long <- convolve_with_hrf(X, hrf_long)
  expect_equal(dim(result_long), dim(X))
  
  # Multiple HRFs (for future CBD support)
  hrf_matrix <- matrix(rnorm(3 * 16), nrow = 3, ncol = 16)
  result_multi <- convolve_with_hrf(X, hrf_matrix)
  expect_equal(dim(result_multi), dim(X))
})

test_that("Gradient computation is numerically stable", {
  skip_if_not(exists("compute_gradient_fista_rcpp"), "Rcpp functions not compiled")
  
  # Test with extreme values
  W <- matrix(c(1e6, 1e-6, 1, 1e6, 1e-6, 1), nrow = 3, ncol = 2)
  X <- matrix(c(1e-6, 1e6), nrow = 2, ncol = 1)
  hrf <- c(1, 0.5, 0.25)
  
  # Compute convolved X
  H_star_X <- stance:::convolve_rows_rcpp(X, hrf)
  
  # This should not produce NaN or Inf
  grad <- stance:::compute_gradient_fista_rcpp(
    Y_or_WtY = t(W) %*% W %*% H_star_X,
    W = W,
    H_star_X = H_star_X,
    hrf_kernel = hrf,
    precomputed_WtY = FALSE
  )
  
  expect_true(all(is.finite(grad)))
})

test_that("Convergence checking works correctly", {
  # Test convergence detection
  values_converging <- c(100, 50, 30, 25, 24, 23.5, 23.4, 23.35, 23.34, 23.34)
  expect_true(stance:::check_convergence(values_converging, tol = 1e-3))
  
  # Test non-convergence
  values_diverging <- c(100, 90, 80, 70, 60, 50, 40, 30, 20, 10)
  expect_false(stance:::check_convergence(values_diverging, tol = 1e-3))
  
  # Test with few values
  values_few <- c(100, 90)
  expect_false(stance:::check_convergence(values_few))
  
  # Test with constant values (should converge)
  values_const <- rep(50, 10)
  expect_true(stance:::check_convergence(values_const))
})

test_that("Total variation computation is correct", {
  skip_if_not(exists("compute_tv_rcpp"), "Rcpp functions not compiled")
  
  # Test simple cases
  X1 <- matrix(c(1, 2, 3, 4, 5), nrow = 1)
  tv1 <- stance:::compute_tv_rcpp(X1)
  expect_equal(tv1, 4)  # |2-1| + |3-2| + |4-3| + |5-4| = 4
  
  # Test with multiple rows
  X2 <- matrix(c(1, 2, 4, 7,
                 3, 5, 8, 10), nrow = 2, byrow = TRUE)
  tv2 <- stance:::compute_tv_rcpp(X2)
  expect_equal(tv2, 6 + 5)  # Row 1: 1+2+3=6, Row 2: 2+3=5
  
  # Test edge cases
  X_const <- matrix(rep(5, 10), nrow = 1)
  expect_equal(stance:::compute_tv_rcpp(X_const), 0)
  
  X_single <- matrix(5, nrow = 1, ncol = 1)
  expect_equal(stance:::compute_tv_rcpp(X_single), 0)
})

test_that("CLD handles extreme parameter values", {
  sim <- simulate_fmri_data(V = 100, T = 50, K = 3)
  
  # Extreme lambda_tv
  cld_high_lambda <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    lambda_tv = 1000,  # Very high regularization
    verbose = FALSE
  )
  expect_s3_class(cld_high_lambda, "ContinuousLinearDecoder")
  
  # Zero lambda_tv
  cld_zero_lambda <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    lambda_tv = 0,  # No regularization
    verbose = FALSE
  )
  expect_s3_class(cld_zero_lambda, "ContinuousLinearDecoder")
  
  # Very small rank
  cld_small_rank <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    rank = 1,
    verbose = FALSE
  )
  expect_s3_class(cld_small_rank, "ContinuousLinearDecoder")
})
