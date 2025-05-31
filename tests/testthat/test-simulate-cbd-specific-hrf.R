library(testthat)

context("simulate_cbd_data with voxel HRFs")

test_that("simulate_cbd_data generates voxel-specific HRFs", {
  basis <- create_hrf_basis_canonical(TR = 0.8, duration_secs = 16, n_derivatives = 1)
  V <- 10
  L <- ncol(basis)
  true_H <- matrix(runif(V * L, -0.2, 0.2), V, L)
  sim <- simulate_cbd_data(V = V, T = 40, K = 2,
                           hrf_basis = basis, true_H = true_H,
                           TR = 0.8, snr = 5, verbose = FALSE)
  expect_equal(dim(sim$H), c(V, L))
  expect_equal(dim(sim$Y), c(V, 40))
  expect_true(all.equal(sim$hrf_basis, basis))
})

test_that("neuroimaging output works with voxel HRFs", {
  basis <- create_hrf_basis_canonical(TR = 1, duration_secs = 16, n_derivatives = 0)
  V <- 27
  L <- ncol(basis)
  true_H <- matrix(rnorm(V * L), V, L)
  sim <- simulate_cbd_data(V = V, T = 20, K = 2,
                           hrf_basis = basis, true_H = true_H,
                           return_neuroim = TRUE, dims = c(3,3,3),
                           verbose = FALSE)
  expect_s4_class(sim$Y_neurovol, "NeuroVec")
  expect_equal(dim(sim$H), c(V, L))
})
