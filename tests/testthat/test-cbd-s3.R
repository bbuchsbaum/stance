library(testthat)
library(stance)

# S3 method tests for ContinuousBayesianDecoder

test_that("print and summary methods work", {
  sim <- simulate_fmri_data(V = 20, T = 30, K = 3, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 3, r = 5)

  pr <- capture.output(print(cbd))
  expect_true(any(grepl("Continuous Bayesian Decoder", pr)))

  sm <- capture.output(summary(cbd))
  expect_true(any(grepl("Average state probabilities", sm)))
})

test_that("coef and fitted methods return expected structures", {
  sim <- simulate_fmri_data(V = 15, T = 25, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 4)

  co <- coef(cbd)
  expect_true(all(c("U", "V", "Pi", "pi0", "sigma2", "H_v") %in% names(co)))

  Yhat <- fitted(cbd)
  expect_true(is.matrix(Yhat))
  expect_equal(dim(Yhat), c(15, 25))
})
