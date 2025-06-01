library(testthat)
library(stance)

# Ensure e_step behaves the same under both engines

test_that("e_step produces consistent posteriors across engines", {
  set.seed(101)
  sim <- simulate_fmri_data(V = 6, T = 10, K = 2, algorithm = "CBD", verbose = FALSE)

  cbd_cpp <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2, engine = "cpp")
  cbd_cpp$.__enclos_env__$private$.e_step()
  gamma_cpp <- cbd_cpp$.__enclos_env__$private$.S_gamma
  xi_cpp <- cbd_cpp$.__enclos_env__$private$.S_xi

  cbd_r <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2, engine = "R")
  cbd_r$.__enclos_env__$private$.e_step()
  gamma_r <- cbd_r$.__enclos_env__$private$.S_gamma
  xi_r <- cbd_r$.__enclos_env__$private$.S_xi

  expect_equal(gamma_cpp, gamma_r, tolerance = 1e-6)
  expect_equal(xi_cpp, xi_r, tolerance = 1e-6)
})
