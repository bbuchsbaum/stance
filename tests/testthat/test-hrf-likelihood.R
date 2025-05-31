library(testthat)
library(stance)

context("Likelihood uses voxel HRFs")

test_that("compute_log_likelihoods reacts to HRF changes", {
  sim <- simulate_cbd_data(V = 5, T = 10, K = 2, verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2,
                                       hrf_basis = sim$hrf_basis, engine = "R")
  ll1 <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
  cbd$.__enclos_env__$private$.H_v <- cbd$.__enclos_env__$private$.H_v * 2
  ll2 <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
  expect_true(max(abs(ll1 - ll2)) > 0)
})
