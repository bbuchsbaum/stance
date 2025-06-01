context("voxel-specific convolution")

library(stance)

test_that("convolve_voxel_hrf_rcpp matches R implementation", {
  set.seed(1)
  K <- 2
  T_len <- 20
  V <- 3
  L <- 4
  design <- matrix(rbinom(K * T_len, 1, 0.3), K, T_len)
  hrfs <- matrix(runif(V * L), V, L)

  r_res <- array(0, dim = c(V, K, T_len))
  for (v in seq_len(V)) {
    conv <- convolve_with_hrf(design, hrfs[v, ])
    r_res[v, , ] <- conv
  }

  cpp_res <- convolve_voxel_hrf_rcpp(design, hrfs, fft_threshold = 1000L, n_threads = 1L)

  expect_equal(cpp_res, r_res, tolerance = 1e-8)
})

test_that("FFT convolution matches direct version", {
  set.seed(2)
  K <- 2; T_len <- 300; V <- 2; L <- 5
  design <- matrix(rnorm(K * T_len), K, T_len)
  hrfs <- matrix(runif(V * L), V, L)

  direct <- convolve_voxel_hrf_rcpp(design, hrfs, fft_threshold = 1000L, n_threads = 1L)
  fft_res <- convolve_voxel_hrf_fft_rcpp(design, hrfs, n_threads = 1L)

  expect_equal(fft_res, direct, tolerance = 1e-8)
})
