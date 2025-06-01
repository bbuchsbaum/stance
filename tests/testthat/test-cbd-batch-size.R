library(testthat)
library(stance)

test_that("fit() errors when batch_size exceeds timepoints", {
  sim <- simulate_fmri_data(V = 10, T = 20, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)
  expect_error(
    cbd$fit(max_iter = 1, verbose = FALSE, batch_size = 25),
    "batch_size must be less than or equal to the number of timepoints"
  )
})
