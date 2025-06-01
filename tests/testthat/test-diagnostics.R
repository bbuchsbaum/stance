library(testthat)
library(stance)

set.seed(1)
sim <- simulate_fmri_data(V = 10, T = 20, K = 2, algorithm = "CBD", verbose = FALSE)
cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 4)
cbd$fit(max_iter = 2, verbose = FALSE)


diag <- cbd$diagnose_model_fit()

test_that("diagnostics list has expected components", {
  expect_true(is.list(diag))
  expect_true("convergence" %in% names(diag))
  expect_true("performance" %in% names(diag))
})
