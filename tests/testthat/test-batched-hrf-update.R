library(testthat)
library(stance)

context("Batched HRF coefficient update")

test_that("update_hrf_coefficients_batched_cpp returns correct dimensions", {
  basis <- create_hrf_basis_canonical(TR = 1, duration_secs = 8, n_derivatives = 0)
  Y <- matrix(rnorm(5 * nrow(basis)), 5, nrow(basis))
  H <- stance:::update_hrf_coefficients_batched_cpp(Y, basis)
  expect_equal(dim(H), c(5, ncol(basis)))
})

test_that("update_hrf_coefficients_batched_cpp handles non-SPD systems", {
  basis <- matrix(1, nrow = 6, ncol = 2)  # singular design
  Y <- matrix(rnorm(4 * 6), 4, 6)
  H <- stance:::update_hrf_coefficients_batched_cpp(Y, basis, ridge = 1e-3)
  expect_equal(dim(H), c(4, 2))
})
