library(testthat)
library(stance)

# Validate that R fallback functions match compiled versions when available

test_that("forward_backward fallback matches C++", {
  set.seed(1)
  Y <- matrix(rnorm(10), 2, 5)
  U <- matrix(rnorm(4), 2, 2)
  V <- matrix(rnorm(4), 2, 2)
  H_v <- matrix(rnorm(2), 2, 1)
  hrf_basis <- matrix(1, 1, 1)
  Pi <- matrix(c(0.7,0.3,0.4,0.6),2,2,byrow=TRUE)
  pi0 <- c(0.6,0.4)
  sigma2 <- 1

  res_cpp <- forward_backward_cpp(Y, U, V, H_v, hrf_basis, Pi, pi0, sigma2)
  res_r   <- forward_backward_r(Y, U, V, H_v, hrf_basis, Pi, pi0, sigma2)

  expect_equal(res_cpp$gamma, res_r$gamma, tolerance = 1e-6)
  expect_equal(res_cpp$xi, res_r$xi, tolerance = 1e-6)
})


test_that("spatial and HRF updates fallback match C++", {
  set.seed(2)
  Y <- matrix(rnorm(12), 3, 4)
  S_gamma <- matrix(runif(6), 2, 3)
  S_gamma <- sweep(S_gamma, 2, colSums(S_gamma), "/")
  U <- matrix(rnorm(9), 3, 3)
  V <- matrix(rnorm(6), 2, 3)
  H_v <- matrix(rnorm(3), 3, 1)
  hrf_basis <- matrix(1, 1, 1)

  sp_cpp <- update_spatial_components_cpp(Y, S_gamma, H_v, hrf_basis, U, V)
  sp_r   <- update_spatial_components_r(Y, S_gamma, H_v, hrf_basis, U, V)

  expect_equal(sp_cpp$U, sp_r$U, tolerance = 1e-6)
  expect_equal(sp_cpp$V, sp_r$V, tolerance = 1e-6)

  hr_cpp <- update_hrf_coefficients_gmrf_cpp(Y, S_gamma, U, V, hrf_basis,
                                             NULL, 1, 1)
  hr_r   <- update_hrf_coefficients_r(Y, S_gamma, U, V, hrf_basis,
                                      NULL, 1, 1)
  expect_equal(hr_cpp, hr_r, tolerance = 1e-6)
})

