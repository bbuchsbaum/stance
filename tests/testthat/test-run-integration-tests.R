library(testthat)
library(stance)

# Test run_integration_tests wrapper

test_that("run_integration_tests executes configurations", {
  sim <- simulate_fmri_data(V = 8, T = 10, K = 2, algorithm = "CBD", verbose = FALSE)
  configs <- list(
    base = list(use_gmrf = FALSE, use_rcpp = FALSE),
    base2 = list(use_gmrf = FALSE, use_rcpp = FALSE)
  )
  res <- run_integration_tests(configs, sim, max_iter = 1)
  expect_type(res, "list")
  expect_equal(names(res), names(configs))
  for (r in res) {
    expect_true(is.numeric(r$time))
    expect_true(is.logical(r$converged))
  }
})
