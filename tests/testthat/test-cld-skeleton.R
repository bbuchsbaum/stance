library(stance)

test_that("ContinuousLinearDecoder class exists", {
  expect_true(exists("ContinuousLinearDecoder"))
  expect_s3_class(ContinuousLinearDecoder, "R6ClassGenerator")
})

test_that("CLD initialization works", {
  # CLD requires Y and S_design parameters
  Y <- matrix(rnorm(100), 10, 10)
  S_design <- matrix(rnorm(30), 3, 10)
  
  cld <- ContinuousLinearDecoder$new(
    Y = Y,
    S_design = S_design,
    verbose = FALSE
  )
  expect_s3_class(cld, "ContinuousLinearDecoder")
})

test_that("CLD has expected public methods", {
  cld_gen <- ContinuousLinearDecoder
  
  # Check public methods exist
  expect_true("initialize" %in% names(cld_gen$public_methods))
  expect_true("fit" %in% names(cld_gen$public_methods))
  expect_true("get_spatial_maps" %in% names(cld_gen$public_methods))
  expect_true("get_activations" %in% names(cld_gen$public_methods))
  expect_true("get_diagnostics" %in% names(cld_gen$public_methods))
  expect_true("print" %in% names(cld_gen$public_methods))
})

test_that("CLD has expected active bindings", {
  cld_gen <- ContinuousLinearDecoder
  
  # Check active bindings exist
  expect_true("W" %in% names(cld_gen$active))
  expect_true("X_hat" %in% names(cld_gen$active))
  expect_true("hrf" %in% names(cld_gen$active))
  expect_true("lambda_tv" %in% names(cld_gen$active))
  expect_true("fitted" %in% names(cld_gen$active))
})

test_that("CLD math helpers work", {
  # Test soft thresholding
  x <- c(-2, -1, 0, 1, 2)
  result <- stance:::soft_threshold(x, lambda = 0.5)
  expect_equal(result, c(-1.5, -0.5, 0, 0.5, 1.5))
  
  # Test compute_tv
  X <- matrix(c(1, 2, 4, 7, 3, 5, 8, 10), nrow = 2, byrow = TRUE)
  tv_val <- stance:::compute_tv(X)
  expect_equal(tv_val, sum(abs(c(1, 2, 3))) + sum(abs(c(2, 3, 2))))
  
  # Test initialize_states
  states_zeros <- stance:::initialize_states(3, 10, "zeros")
  expect_equal(dim(states_zeros), c(3, 10))
  expect_true(all(states_zeros == 0))
  
  states_uniform <- stance:::initialize_states(3, 10, "uniform")
  expect_true(all(states_uniform == 1/3))
})

test_that("check_convergence works", {
  # Not converged
  values1 <- c(100, 90, 80, 70, 60)
  expect_false(stance:::check_convergence(values1))
  
  # Converged (very small relative change)
  values2 <- c(100, 50, 25, 24.9, 24.85, 24.845, 24.844)
  # The relative change between 24.9 and 24.844 is (24.9-24.844)/24.9 = 0.00225 > 1e-4
  # So this should not converge with default tolerance
  expect_false(stance:::check_convergence(values2))
  
  # Actually converged values
  values3 <- c(100, 50, 25.0, 24.9999, 24.9998, 24.99975, 24.9997)
  # Relative change is (25.0 - 24.9997)/25.0 = 0.000012 < 1e-4
  expect_true(stance:::check_convergence(values3))
  
  # Too few values
  values4 <- c(100, 90)
  expect_false(stance:::check_convergence(values4))
})
