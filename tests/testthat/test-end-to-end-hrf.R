library(testthat)
library(stance)

context("CBD end-to-end with voxel HRFs")

set.seed(123)

basis <- create_hrf_basis_canonical(TR = 0.8, duration_secs = 16, n_derivatives = 1)
V <- 20
T_len <- 50
K <- 2
L <- ncol(basis)
true_H <- matrix(rnorm(V * L, sd = 0.1), V, L)

sim <- simulate_cbd_data(V = V, T = T_len, K = K, hrf_basis = basis,
                         true_H = true_H, TR = 0.8, snr = 3, verbose = FALSE)

cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = K, r = 3,
                                     hrf_basis = basis, engine = "cpp")

# Fit a few VB iterations
fit_time <- system.time({
  cbd$fit(max_iter = 5, verbose = FALSE)
})

est_H <- cbd$get_hrf_estimates()
cor_vals <- numeric(5)
for (i in seq_len(5)) {
  cor_vals[i] <- cor(est_H[i, ], true_H[i, ])
}

conv <- cbd$get_convergence()

# speed benchmark comparing R and C++ likelihood kernels
Y_proj <- cbd$.__enclos_env__$private$.Y_proj
Vmat <- cbd$.__enclos_env__$private$.V
hrf_kernel <- as.vector(basis %*% colMeans(est_H))
noise_var <- cbd$.__enclos_env__$private$.sigma2

r_time <- system.time({
  stance:::compute_log_likelihoods_r(Y_proj, Vmat, hrf_kernel, noise_var)
})
cpp_time <- system.time({
  compute_log_likelihoods_rcpp(Y_proj, Vmat, hrf_kernel, noise_var)
})

test_that("HRF recovery correlation is reasonable", {
  expect_true(mean(abs(cor_vals)) > 0.3)
})

test_that("ELBO increases during fitting", {
  expect_true(all(diff(conv$elbo_history) >= -1e-8))
})

test_that("Rcpp likelihood is faster than R version", {
  expect_true(r_time["elapsed"] / cpp_time["elapsed"] >= 5)
})

# Neuroim2 compatibility

test_that("pipeline works with NeuroVec input", {
  skip_if_not_installed("neuroim2")
  sim_ni <- simulate_cbd_data(V = V, T = T_len, K = K, hrf_basis = basis,
                              true_H = true_H, TR = 0.8, snr = 3,
                              return_neuroim = TRUE, dims = c(4,5,1),
                              verbose = FALSE)
  cbd_ni <- ContinuousBayesianDecoder$new(Y = sim_ni$Y_neurovol,
                                          K = K, r = 3,
                                          hrf_basis = basis, engine = "cpp")
  cbd_ni$fit(max_iter = 3, verbose = FALSE)
  expect_equal(dim(cbd_ni$get_hrf_estimates()), c(V, L))
})

