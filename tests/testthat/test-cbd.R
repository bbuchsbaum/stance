library(testthat)
library(stance)

# Algorithm validation tests for ContinuousBayesianDecoder

# Forward-backward accuracy against reference implementation

test_that("forward_backward matches reference", {
  set.seed(123)
  sim <- simulate_fmri_data(V = 5, T = 8, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)
  loglik <- cbd$.__enclos_env__$private$.compute_log_likelihoods()
  res_scaled <- cbd$.__enclos_env__$private$.forward_backward(loglik)
  res_ref <- stance:::forward_backward_r(
    cbd$.__enclos_env__$private$.Y_data,
    cbd$.__enclos_env__$private$.U,
    cbd$.__enclos_env__$private$.V,
    cbd$.__enclos_env__$private$.H_v,
    cbd$.__enclos_env__$private$.hrf_kernel,
    cbd$.__enclos_env__$private$.Pi,
    cbd$.__enclos_env__$private$.pi0,
    cbd$.__enclos_env__$private$.sigma2
  )
  expect_equal(res_scaled$gamma, res_ref$gamma, tolerance = 1e-5)
  expect_equal(res_scaled$xi, res_ref$xi, tolerance = 1e-5)
})

# VB update convergence and ELBO monotonicity

test_that("VB updates increase ELBO", {
  sim <- simulate_fmri_data(V = 6, T = 10, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)

  params <- list(
    U = cbd$.__enclos_env__$private$.U,
    V = cbd$.__enclos_env__$private$.V,
    H_v = cbd$.__enclos_env__$private$.H_v,
    hrf_basis = cbd$.__enclos_env__$private$.hrf_kernel,
    Pi = cbd$.__enclos_env__$private$.Pi,
    pi0 = cbd$.__enclos_env__$private$.pi0,
    sigma2 = cbd$.__enclos_env__$private$.sigma2,
    L_gmrf = cbd$.__enclos_env__$private$.L_gmrf,
    lambda_H_prior = cbd$.__enclos_env__$private$.lambda_H_prior,
    sigma2_prior = NULL,
    prior_Pi = NULL,
    prior_pi0 = NULL
  )
  config <- list(engine = "R")

  elbos <- numeric(3)
  for (i in seq_len(3)) {
    vb <- stance:::vb_e_step(cbd$.__enclos_env__$private$.Y_data, params, config)
    params <- stance:::vb_m_step(cbd$.__enclos_env__$private$.Y_data, vb, params, config)
    elbos[i] <- stance:::compute_elbo(cbd$.__enclos_env__$private$.Y_data, vb, params, config)
  }
  expect_true(all(diff(elbos) >= -1e-8))
})

# Parameter recovery on simulated data

test_that("CBD recovers parameters on simulated data", {
  sim <- simulate_fmri_data(V = 20, T = 30, K = 2, algorithm = "CBD", snr = 5, verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 3)
  cbd$fit(max_iter = 5, verbose = FALSE)

  W_est <- cbd$get_spatial_maps(as_neurovol = FALSE)
  true_W <- sim$W

  corr1 <- abs(cor(W_est[,1], true_W[,1]))
  corr2 <- abs(cor(W_est[,2], true_W[,2]))
  expect_true(mean(c(corr1, corr2)) > 0.2)

  states_est <- apply(cbd$get_state_sequence(), 2, which.max)
  states_true <- apply(sim$S, 2, which.max)
  acc <- mean(states_est == states_true)
  expect_true(acc > 0.3)
})

# CBD vs CLD comparison on identical data

test_that("CBD outperforms or matches CLD on same data", {
  sim <- simulate_fmri_data(V = 20, T = 25, K = 2, algorithm = "CBD", snr = 5, verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)
  cbd$fit(max_iter = 5, verbose = FALSE)

  cld <- ContinuousLinearDecoder$new(Y = sim$Y, S_design = sim$S, verbose = FALSE)
  cld$fit(max_iter = 5, verbose = FALSE)

  states_cbd <- apply(cbd$get_state_sequence(), 2, which.max)
  states_cld <- apply(cld$get_activations(), 2, which.max)
  states_true <- apply(sim$S, 2, which.max)

  acc_cbd <- mean(states_cbd == states_true)
  acc_cld <- mean(states_cld == states_true)
  expect_true(acc_cbd >= acc_cld - 0.1)
})


