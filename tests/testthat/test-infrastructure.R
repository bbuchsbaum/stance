# Test shared infrastructure components

test_that("package loads correctly", {
  expect_true(requireNamespace("stance", quietly = TRUE))
})

test_that("shared utilities are available", {
  # Data structure utilities
  expect_true(exists("validate_fmri_input"))
  expect_true(exists("extract_data_matrix"))
  expect_true(exists("restore_spatial_structure"))
  expect_true(exists("check_temporal_alignment"))
  
  # HRF utilities
  expect_true(exists("setup_hrf_kernel"))
  expect_true(exists("convolve_with_hrf"))
  expect_true(exists("hrf_basis_matrix"))
  expect_true(exists("validate_hrf_spec"))
  
  # Simulation utilities
  expect_true(exists("simulate_fmri_data"))
  expect_true(exists("generate_event_design"))
  
  # Visualization utilities
  expect_true(exists("plot_spatial_maps"))
  expect_true(exists("plot_state_timecourse"))
  expect_true(exists("plot_convergence"))
})

test_that("data structure conversion preserves information", {
  # Create test data
  V <- 100
  T <- 50
  Y_mat <- matrix(rnorm(V * T), V, T)
  
  # Matrix -> validated -> matrix roundtrip
  validated <- validate_fmri_input(Y_mat)
  expect_equal(validated$data, Y_mat)
  expect_equal(validated$dims["V"], V)
  expect_equal(validated$dims["T"], T)
  
  # Extract and restore
  extracted <- extract_data_matrix(Y_mat)
  expect_equal(extracted$data, Y_mat)
})

test_that("HRF utilities integrate with fmrireg", {
  # Test HRF setup
  hrf1 <- setup_hrf_kernel("spmg1", TR = 2)
  expect_true(is.numeric(hrf1))
  expect_true(length(hrf1) > 1)
  
  # Test normalization
  expect_equal(sum(hrf1^2), 1, tolerance = 1e-10)
  
  # Test basis matrix generation
  basis <- hrf_basis_matrix("spmg2", TR = 2, len = 20)
  expect_true(is.matrix(basis))
  expect_equal(ncol(basis), 2)  # Canonical + temporal derivative
})

test_that("simulation framework supports both algorithms", {
  # CLD simulation
  sim_cld <- simulate_fmri_data(V = 100, T = 50, K = 2, algorithm = "CLD")
  expect_equal(sim_cld$params$algorithm, "CLD")
  expect_true(any(sim_cld$S > 0 & sim_cld$S < 1))  # Continuous values
  
  # CBD simulation
  sim_cbd <- simulate_fmri_data(V = 100, T = 50, K = 2, algorithm = "CBD")
  expect_equal(sim_cbd$params$algorithm, "CBD")
  expect_true(all(sim_cbd$S %in% c(0, 1)))  # Binary values
  
  # Both should have same structure
  expect_equal(names(sim_cld), names(sim_cbd))
})

test_that("visualization functions handle various inputs", {
  # Create test data
  W <- matrix(rnorm(100 * 3), 100, 3)
  X <- matrix(rnorm(3 * 50), 3, 50)
  
  # These should not error
  expect_silent({
    pdf(NULL)  # Null device
    plot_spatial_maps(W)
    plot_state_timecourse(X)
    plot_convergence(c(100, 90, 80, 75, 74, 73))
    dev.off()
  })
})

test_that("integration between components works", {
  # Full pipeline test
  # 1. Simulate data
  sim <- simulate_fmri_data(V = 64, T = 100, K = 2, algorithm = "CLD")
  
  # 2. Validate input
  validated <- validate_fmri_input(sim$Y)
  expect_equal(dim(validated$data), dim(sim$Y))
  
  # 3. HRF processing
  hrf <- setup_hrf_kernel("spmg1")
  X_conv <- convolve_with_hrf(sim$S, hrf)
  expect_equal(dim(X_conv), dim(sim$S))
  
  # 4. Check alignment
  design <- matrix(1, nrow = 2, ncol = 100)
  expect_true(check_temporal_alignment(validated, design, TR = 2))
})

test_that("edge cases are handled gracefully", {
  # Empty data
  expect_error(validate_fmri_input(matrix(nrow = 0, ncol = 0)))
  
  # Single voxel
  Y_single <- matrix(rnorm(50), 1, 50)
  val_single <- validate_fmri_input(Y_single)
  expect_equal(val_single$dims["V"], 1)
  
  # Single timepoint
  Y_time <- matrix(rnorm(100), 100, 1)
  val_time <- validate_fmri_input(Y_time)
  expect_equal(val_time$dims["T"], 1)
  
  # Invalid HRF
  expect_error(setup_hrf_kernel("invalid_hrf"))
  expect_error(validate_hrf_spec(NULL))
})

test_that("performance characteristics are reasonable", {
  # Time key operations
  V <- 1000
  T <- 200
  K <- 3
  
  # Data validation should be fast
  Y <- matrix(rnorm(V * T), V, T)
  time_validate <- system.time({
    validated <- validate_fmri_input(Y)
  })
  expect_true(time_validate["elapsed"] < 0.1)  # Should be < 100ms
  
  # HRF convolution should scale linearly
  X <- matrix(rnorm(K * T), K, T)
  hrf <- setup_hrf_kernel("spmg1")
  time_conv <- system.time({
    X_conv <- convolve_with_hrf(X, hrf)
  })
  expect_true(time_conv["elapsed"] < 0.5)  # Should be < 500ms
})