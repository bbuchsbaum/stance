library(testthat)
library(stance)

# Test log likelihood calculation

test_that("compute_log_likelihoods returns higher likelihood for true state", {
  set.seed(1)
  sim <- simulate_fmri_data(V = 10, T = 20, K = 3, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 3, r = 4)
  loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()

  expect_equal(dim(loglik), c(3, 20))
  true_states <- apply(sim$S, 2, function(x) which(x == 1))

  higher <- vapply(seq_len(20), function(t) {
    ll_t <- loglik[, t]
    ll_t[true_states[t]] >= max(ll_t)
  }, logical(1))

  expect_true(all(higher))
})
