library(stance)

test_that("setup_hrf_kernel handles numeric input", {
  # Test with numeric vector
  hrf_numeric <- c(0, 0.1, 0.5, 0.8, 0.6, 0.3, 0.1, 0)
  
  # Without normalization
  result1 <- setup_hrf_kernel(hrf_numeric, normalize = FALSE)
  expect_equal(result1, hrf_numeric)
  
  # With normalization
  result2 <- setup_hrf_kernel(hrf_numeric, normalize = TRUE)
  expect_equal(sum(result2^2), 1, tolerance = 1e-10)
  
  # Check proportions are maintained
  expect_equal(result2[3] / result2[4], hrf_numeric[3] / hrf_numeric[4])
})

test_that("setup_hrf_kernel handles character specifications", {
  # Test canonical HRF
  hrf1 <- setup_hrf_kernel("spmg1", TR = 2, len = 32)
  expect_true(is.numeric(hrf1))
  expect_true(length(hrf1) > 0)
  expect_equal(sum(hrf1^2), 1, tolerance = 1e-10)
  
  # Test aliases
  hrf_canonical <- setup_hrf_kernel("canonical", TR = 2, len = 32)
  hrf_spm <- setup_hrf_kernel("spm", TR = 2, len = 32)
  expect_equal(hrf_canonical, hrf_spm)
})

test_that("validate_hrf_spec catches invalid inputs", {
  # NULL input
  expect_error(validate_hrf_spec(NULL), "cannot be NULL")
  
  # Invalid character
  expect_error(validate_hrf_spec("invalid_hrf"), "not a recognized HRF type")
  
  # Too short numeric
  expect_error(validate_hrf_spec(c(1)), "must have length >= 2")
  
  # Non-finite values
  expect_error(validate_hrf_spec(c(1, NA, 3)), "contains non-finite values")
  
  # All zeros
  expect_error(validate_hrf_spec(c(0, 0, 0)), "is all zeros")
  
  # Valid inputs should pass
  expect_true(validate_hrf_spec("spmg1"))
  expect_true(validate_hrf_spec(c(0.1, 0.5, 0.3)))
})

test_that("convolve_with_hrf works correctly", {
  # Simple test signal
  X <- matrix(0, nrow = 2, ncol = 10)
  X[1, 3] <- 1  # Impulse at t=3
  X[2, 5] <- 1  # Impulse at t=5
  
  # Simple HRF
  hrf <- c(0.2, 0.5, 0.3)
  
  # Direct convolution
  result <- convolve_with_hrf(X, hrf, method = "direct")
  
  # Check dimensions
  expect_equal(dim(result), dim(X))
  
  # Check convolution results
  expect_equal(result[1, 3:5], hrf)
  expect_equal(result[2, 5:7], hrf)
  
  # Test with FFT
  result_fft <- convolve_with_hrf(X, hrf, method = "fft")
  expect_equal(result, result_fft, tolerance = 1e-10)
})

test_that("convolve_with_hrf handles automatic method selection", {
  # Small signal - should use direct
  X_small <- matrix(rnorm(3 * 100), nrow = 3, ncol = 100)
  hrf <- setup_hrf_kernel("spmg1", TR = 1, len = 20)
  
  result1 <- convolve_with_hrf(X_small, hrf, use_fft_threshold = 256)
  
  # Large signal - should use FFT
  X_large <- matrix(rnorm(3 * 500), nrow = 3, ncol = 500)
  result2 <- convolve_with_hrf(X_large, hrf, use_fft_threshold = 256)
  
  # Both should produce valid results
  expect_equal(dim(result1), dim(X_small))
  expect_equal(dim(result2), dim(X_large))
})

test_that("hrf_basis_matrix generates correct bases", {
  # Test SPM bases
  basis1 <- hrf_basis_matrix("spmg1", TR = 2, len = 32)
  expect_true(is.matrix(basis1))
  expect_equal(ncol(basis1), 1)  # Single basis
  
  basis2 <- hrf_basis_matrix("spmg2", TR = 2, len = 32)
  expect_equal(ncol(basis2), 2)  # Canonical + temporal derivative
  
  basis3 <- hrf_basis_matrix("spmg3", TR = 2, len = 32)
  expect_equal(ncol(basis3), 3)  # Canonical + temporal + dispersion
  
  # Check attributes
  expect_equal(attr(basis3, "basis_type"), "spmg3")
  expect_equal(attr(basis3, "TR"), 2)
})

test_that("hrf_convolution_matrix creates correct Toeplitz structure", {
  hrf <- c(0.3, 0.5, 0.2)
  n_time <- 10
  
  # Sparse matrix
  conv_mat_sparse <- hrf_convolution_matrix(hrf, n_time, sparse = TRUE)
  expect_s4_class(conv_mat_sparse, "sparseMatrix")
  expect_equal(dim(conv_mat_sparse), c(n_time, n_time))
  
  # Dense matrix
  conv_mat_dense <- hrf_convolution_matrix(hrf, n_time, sparse = FALSE)
  expect_true(is.matrix(conv_mat_dense))
  expect_equal(dim(conv_mat_dense), c(n_time, n_time))

  # Sparse and dense representations should be identical element-wise
  expect_equal(as.matrix(conv_mat_sparse), conv_mat_dense)
  
  # Check Toeplitz structure
  # First column should have HRF at the beginning
  expect_equal(as.numeric(conv_mat_dense[1:3, 1]), hrf)
  
  # Check that convolution with delta gives HRF
  delta <- rep(0, n_time)
  delta[1] <- 1
  result <- as.numeric(conv_mat_dense %*% delta)
  expect_equal(result[1:3], hrf)
})

test_that("convolve_with_hrf_transposed performs transposed convolution", {
  # Create test signal
  X <- matrix(rnorm(2 * 20), nrow = 2, ncol = 20)
  hrf <- c(0.2, 0.5, 0.3)
  
  # Transposed convolution
  result <- convolve_with_hrf_transposed(X, hrf)
  
  # Should have same dimensions
  expect_equal(dim(result), dim(X))
  
  # Compare with manual computation using reversed HRF
  expected <- convolve_with_hrf(X, rev(hrf))
  expect_equal(result, expected)
})

test_that("convolve_fft internal function works correctly", {
  # Simple vectors
  x <- c(1, 0, 0, 0, 0)
  h <- c(0.5, 0.3, 0.2)
  
  # This function is internal but we can test it
  result <- stance:::convolve_fft(x, h)
  
  expect_equal(length(result), length(x))
  expect_equal(result[1:3], h, tolerance = 1e-10)
})

test_that("create_hrf_basis_canonical generates orthonormal basis", {
  basis <- create_hrf_basis_canonical()
  expect_true(is.matrix(basis))
  expect_equal(ncol(basis), 3)
  expect_equal(nrow(basis), floor(32 / 0.8) + 1)
  cross <- crossprod(basis)
  expect_equal(cross, diag(ncol(basis)), tolerance = 1e-6)
})

test_that("create_hrf_basis_fir handles TR drift edge cases", {
  basis <- create_hrf_basis_fir(TR = 0.7, duration_secs = 31, fir_resolution_secs = 1)
  expect_true(is.matrix(basis))
  expect_equal(nrow(basis), floor(31 / 0.7) + 1)
  expect_true(ncol(basis) >= 1)
})

test_that("create_hrf_basis_neuroim2 delegates to hrf_basis_matrix", {
  basis_neuro <- create_hrf_basis_neuroim2("spmg3", TR = 1, duration = 16)
  basis_ref   <- hrf_basis_matrix("spmg3", TR = 1, len = 16)
  expect_equal(basis_neuro, basis_ref)
  expect_equal(attr(basis_neuro, "basis_type"), "spmg3")
  expect_equal(attr(basis_neuro, "TR"), 1)
})
