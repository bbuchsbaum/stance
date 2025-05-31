library(testthat)
library(stance)

# Forward-backward algorithm unit test

test_that("forward_backward returns normalized posteriors", {
  set.seed(42)
  sim <- simulate_fmri_data(V = 10, T = 15, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 3)

  loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
  fb <- cbd$.__enclos_env__$private$.forward_backward(loglik)

  gamma <- fb$gamma
  xi <- fb$xi

  expect_equal(dim(gamma), c(2, 15))
  expect_equal(dim(xi), c(2, 2, 14))
  expect_true(all(abs(colSums(gamma) - 1) < 1e-6))
  expect_true(all(abs(apply(xi, 3, sum) - 1) < 1e-6))
  expect_true(is.numeric(fb$log_likelihood))
})
