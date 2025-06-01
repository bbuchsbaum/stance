# Test CLD-specific functionality

test_that("CLD initialization works with simulated data", {
  # Generate simple test data
  sim <- simulate_fmri_data(V = 100, T = 50, K = 3, algorithm = "CLD")
  
  # Create CLD object
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    hrf = "spmg1",
    rank = 10,
    lambda_tv = 0.01,
    verbose = FALSE
  )
  
  expect_s3_class(cld, "ContinuousLinearDecoder")
  expect_equal(cld$lambda_tv, 0.01)
  expect_false(cld$fitted)
})

test_that("CLD handles different input formats", {
  sim <- simulate_fmri_data(V = 64, T = 50, K = 2)
  
  # Matrix input
  cld1 <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    verbose = FALSE
  )
  expect_s3_class(cld1, "ContinuousLinearDecoder")
  
  # Data frame input (should be converted)
  Y_df <- as.data.frame(sim$Y)
  cld2 <- ContinuousLinearDecoder$new(
    Y = Y_df,
    S_design = sim$S,
    verbose = FALSE
  )
  expect_s3_class(cld2, "ContinuousLinearDecoder")
})

test_that("CLD parameter validation works", {
  sim <- simulate_fmri_data(V = 100, T = 50, K = 3)
  
  # Invalid lambda_tv
  expect_error(
    ContinuousLinearDecoder$new(
      Y = sim$Y,
      S_design = sim$S,
      lambda_tv = -1
    ),
    "lambda_tv must be non-negative"
  )
  
  # Mismatched dimensions
  expect_error(
    ContinuousLinearDecoder$new(
      Y = sim$Y,
      S_design = matrix(1, 3, 30)  # Wrong T dimension
    ),
    "S_design has 30 columns but Y has 50 timepoints"
  )
  
  # Rank too high (should warn)
  expect_warning(
    ContinuousLinearDecoder$new(
      Y = sim$Y,
      S_design = sim$S,
      rank = 200  # > min(V, K)
    ),
    "Rank 200 exceeds maximum"
  )
})

test_that("GLM+SVD learning produces reasonable results", {
  # Generate data with known structure
  V <- 200
  T <- 100
  K <- 3
  rank <- 5
  
  # Create low-rank structure
  U_true <- qr.Q(qr(matrix(rnorm(V * rank), V, rank)))
  V_true <- matrix(rnorm(K * rank), K, rank)
  W_true <- U_true %*% t(V_true)
  
  # Generate states
  S_true <- matrix(0, K, T)
  for (k in 1:K) {
    S_true[k, ((k-1)*30 + 1):(k*30)] <- 1
  }
  
  # Generate data
  hrf <- setup_hrf_kernel("spmg1", TR = 1)
  X_true <- convolve_with_hrf(S_true, hrf)
  Y <- W_true %*% X_true + matrix(rnorm(V * T, sd = 0.1), V, T)
  
  # Fit CLD
  cld <- ContinuousLinearDecoder$new(
    Y = Y,
    S_design = S_true,
    hrf = hrf,
    rank = rank,
    verbose = FALSE
  )
  
  # Check learned W
  W_learned <- cld$W
  expect_equal(dim(W_learned), c(V, K))
  
  # Should have reasonable correlation with true W
  # (Note: sign ambiguity in SVD, so check absolute correlation)
  for (k in 1:K) {
    cor_k <- abs(cor(W_true[, k], W_learned[, k]))
    expect_true(cor_k > 0.7)  # Should be reasonably correlated
  }
})

test_that("CLD getters work correctly", {
  sim <- simulate_fmri_data(V = 64, T = 50, K = 2)
  
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    verbose = FALSE
  )
  
  # Get spatial maps
  W <- cld$get_spatial_maps()
  expect_true(is.matrix(W))
  expect_equal(dim(W), c(64, 2))
  
  # Get activations (before fitting)
  X <- cld$get_activations()
  expect_true(is.matrix(X))
  expect_equal(dim(X), c(2, 50))
  expect_true(all(X == 0))  # Should be zeros before fitting
  
  # Get diagnostics
  diag <- cld$get_diagnostics()
  expect_true(is.list(diag))
  expect_true("objective_values" %in% names(diag))
  expect_true("converged" %in% names(diag))
})

test_that("CLD print method works", {
  sim <- simulate_fmri_data(V = 100, T = 50, K = 3)
  
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    verbose = FALSE
  )
  
  # Capture print output
  output <- capture.output(print(cld))
  expect_true(any(grepl("Continuous Linear Decoder", output)))
  expect_true(any(grepl("100 voxels x 50 timepoints", output)))
  expect_true(any(grepl("Number of states: 3", output)))
})

test_that("CLD active bindings work", {
  sim <- simulate_fmri_data(V = 64, T = 50, K = 2)
  
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    lambda_tv = 0.01,
    verbose = FALSE
  )
  
  # Read-only bindings
  expect_true(is.matrix(cld$W))
  expect_true(is.matrix(cld$X_hat))
  expect_true(is.numeric(cld$hrf))
  
  # Read-write binding
  expect_equal(cld$lambda_tv, 0.01)
  cld$lambda_tv <- 0.02
  expect_equal(cld$lambda_tv, 0.02)
  
  # Invalid lambda should error
  expect_error(cld$lambda_tv <- -1, "lambda_tv must be non-negative")
})

test_that("as.list.ContinuousLinearDecoder works", {
  sim <- simulate_fmri_data(V = 64, T = 50, K = 2)
  
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    verbose = FALSE
  )
  
  cld_list <- as.list(cld)
  expect_true(is.list(cld_list))
  expect_true(all(c("W", "X_hat", "hrf", "lambda_tv", "fitted", "diagnostics") %in% names(cld_list)))
})

# Note: FISTA fitting tests would require working Rcpp code compilation
# These are placeholder tests that would work once the package is built

test_that("CLD fitting works (if Rcpp compiled)", {
  skip_if_not(exists("fista_tv_rcpp"), "Rcpp functions not compiled")
  
  sim <- simulate_fmri_data(V = 100, T = 50, K = 2, snr = 5)
  
  cld <- ContinuousLinearDecoder$new(
    Y = sim$Y,
    S_design = sim$S,
    lambda_tv = 0.01,
    verbose = FALSE
  )
  
  # Fit model
  cld$fit(max_iter = 50, tol = 1e-3)
  
  expect_true(cld$fitted)
  expect_true(length(cld$get_diagnostics()$objective_values) > 0)
  
  # Check that states are non-zero after fitting
  X_hat <- cld$get_activations()
  expect_true(any(X_hat != 0))
})

test_that("fista_tv_rcpp results are thread independent", {
  skip_if_not(exists("fista_tv_rcpp"), "Rcpp functions not compiled")

  set.seed(42)
  V <- 10; T_len <- 30; K <- 2
  W <- matrix(rnorm(V * K), V, K)
  X_true <- matrix(rnorm(K * T_len), K, T_len)
  hrf <- c(0.2, 0.5, 0.3)
  Y <- W %*% stance:::convolve_rows_rcpp(X_true, hrf, n_threads = 1L)

  WtY <- t(W) %*% Y
  L <- stance:::estimate_lipschitz_rcpp(W, hrf)

  init <- matrix(0, K, T_len)
  res1 <- stance:::fista_tv_rcpp(WtY, W, hrf, 0.01, L, init, max_iter = 10L, tol = 1e-3, verbose = FALSE, n_threads = 1L)
  res2 <- stance:::fista_tv_rcpp(WtY, W, hrf, 0.01, L, init, max_iter = 10L, tol = 1e-3, verbose = FALSE, n_threads = 2L)

  expect_equal(res1$X_hat, res2$X_hat, tolerance = 1e-8)
})

test_that("CLD handles edge cases gracefully", {
  # Single state
  Y <- matrix(rnorm(100 * 50), 100, 50)
  S <- matrix(1, 1, 50)
  
  cld1 <- ContinuousLinearDecoder$new(Y = Y, S_design = S, verbose = FALSE)
  expect_s3_class(cld1, "ContinuousLinearDecoder")
  
  # Very short time series
  Y2 <- matrix(rnorm(100 * 10), 100, 10)
  S2 <- matrix(rnorm(3 * 10), 3, 10)
  
  cld2 <- ContinuousLinearDecoder$new(Y = Y2, S_design = S2, verbose = FALSE)
  expect_s3_class(cld2, "ContinuousLinearDecoder")
})
