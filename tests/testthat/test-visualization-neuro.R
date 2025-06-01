library(testthat)
library(stance)

context("Enhanced spatial visualization")

test_that("plot_spatial_maps works with decoder object", {
  sim <- simulate_fmri_data(V = 8, T = 10, K = 2, algorithm = "CBD", verbose = FALSE)
  mask <- array(TRUE, c(2,2,2))
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2,
                                       use_gmrf = TRUE, mask = mask)
  pdf(NULL)
  expect_silent(plot_spatial_maps(cbd, states = 1))
  dev.off()
})
