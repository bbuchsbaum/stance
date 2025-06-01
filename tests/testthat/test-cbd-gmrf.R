library(testthat)
library(stance)

test_that("CBD initializes with GMRF prior", {
  mask <- array(TRUE, dim = c(2,2,1))
  Y <- matrix(rnorm(4*10), 4, 10)
  cbd <- ContinuousBayesianDecoder$new(Y = Y, K = 2, r = 2,
                                       use_gmrf = TRUE, lambda_h = 2,
                                       mask = mask)
  expect_true(cbd$.__enclos_env__$private$.use_gmrf)
  L <- cbd$.__enclos_env__$private$.L_gmrf
  expect_s4_class(L, "dgCMatrix")
  expect_equal(dim(L), c(4,4))
})

test_that("GMRF HRF update runs", {
  sim <- simulate_fmri_data(V = 4, T = 6, K = 2, algorithm = "CBD", verbose = FALSE)
  mask <- array(TRUE, dim = c(2,2,1))
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2,
                                       use_gmrf = TRUE, mask = mask)
  cbd$.__enclos_env__$private$.update_hrf_coefficients()
  H <- cbd$get_hrf_estimates()
  expect_equal(dim(H)[1], 4)
})
