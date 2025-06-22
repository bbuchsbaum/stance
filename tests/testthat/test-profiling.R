library(testthat)
library(stance)

# Test profiling of CBD components

test_that("profile_cbd_components returns named timings", {
  sim <- simulate_fmri_data(V = 4, T = 6, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)
  cbd$fit(max_iter = 1, verbose = FALSE)
  timings <- profile_cbd_components(cbd)
  expect_true(is.list(timings))
  expect_named(timings, c("log_lik", "forward_backward", "update_Pi", "update_UV", "update_sigma2"))
  expect_true(all(unlist(timings) >= 0))
})

# Test optimisation suggestions message

test_that("suggest_optimizations highlights slowest components", {
  timings <- list(log_lik = 1, forward_backward = 3, update_Pi = 0.5)
  expect_message(suggest_optimizations(timings, "roi"), "forward_backward")
})

# Test profile_vb_iteration and check_performance_targets

test_that("profile_vb_iteration summarises run and checks targets", {
  skip_if_not_installed("lobstr")
  sim <- simulate_fmri_data(V = 4, T = 6, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)
  cbd$fit(max_iter = 1, verbose = FALSE)
  prof <- profile_vb_iteration(cbd, target = "roi")
  expect_true(all(c("timings", "memory", "meets_target") %in% names(prof)))
  expect_type(prof$meets_target, "logical")
  chk <- check_performance_targets(cbd)
  expect_type(chk, "logical")
})
