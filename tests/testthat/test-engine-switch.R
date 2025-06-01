library(testthat)
library(stance)

# Engine toggling equivalence test

test_that("cpp and R engines produce same forward/backward", {
  set.seed(99)
  sim <- simulate_fmri_data(V = 6, T = 12, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2, engine = "R")

  loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()

  cbd$.__enclos_env__$private$.engine <- "cpp"
  res_cpp <- cbd$.__enclos_env__$private$.forward_backward(loglik)

  cbd$.__enclos_env__$private$.engine <- "R"
  res_r <- cbd$.__enclos_env__$private$.forward_backward(loglik)

  expect_equal(res_cpp$gamma, res_r$gamma, tolerance = 1e-6)
  expect_equal(res_cpp$xi, res_r$xi, tolerance = 1e-6)
})
