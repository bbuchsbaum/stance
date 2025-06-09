library(testthat)
library(stance)

# Test compute_gmrf_prior_term

test_that("compute_gmrf_prior_term matches analytic quadratic form", {
  H_v <- matrix(c(1, 2), nrow = 2, ncol = 1)
  L <- matrix(c(2, -1, -1, 2), nrow = 2, byrow = TRUE)
  lambda <- 2
  result <- stance:::compute_gmrf_prior_term(H_v, L, lambda)
  expected <- -0.5 * lambda * as.numeric(t(H_v[, 1]) %*% L %*% H_v[, 1])
  expect_equal(result, expected)
})

# Test update_hmm_parameters with priors

test_that("update_hmm_parameters applies Dirichlet priors correctly", {
  S_gamma <- matrix(c(0.6, 0.4,
                      0.7, 0.3,
                      0.5, 0.5),
                    nrow = 2, byrow = FALSE)
  S_xi <- array(0, dim = c(2, 2, 2))
  S_xi[, , 1] <- matrix(c(0.7, 0.3,
                          0.4, 0.6),
                        nrow = 2, byrow = TRUE)
  S_xi[, , 2] <- matrix(c(0.5, 0.5,
                          0.2, 0.8),
                        nrow = 2, byrow = TRUE)
  prior_Pi <- matrix(c(2, 1,
                       1, 2),
                     nrow = 2, byrow = TRUE)
  prior_pi0 <- c(2, 1)
  res <- stance:::update_hmm_parameters(S_gamma, S_xi,
                                        prior_Pi = prior_Pi,
                                        prior_pi0 = prior_pi0)

  expected_pi0 <- c(0.8, 0.2)
  expected_Pi <- matrix(c(0.7333333, 0.2666667,
                          0.2,       0.8),
                        nrow = 2, byrow = TRUE)

  expect_equal(res$pi0, expected_pi0, tolerance = 1e-8)
  expect_equal(res$Pi, expected_Pi, tolerance = 1e-8)
  expect_equal(rowSums(res$Pi), rep(1, 2))
})

# Test residual sum of squares and noise variance updates

test_that("noise variance updates use residual sum of squares", {
  U <- matrix(c(1, 1), nrow = 2, ncol = 1)
  V <- matrix(1, nrow = 1, ncol = 1)
  S_gamma <- matrix(c(1, 2, 3), nrow = 1)
  H_v <- matrix(1, nrow = 2, ncol = 1)
  hrf_basis <- matrix(1, nrow = 1, ncol = 1)

  HX <- convolve_with_hrf(S_gamma, 1)
  Y_hat <- U %*% t(V) %*% HX
  Y <- Y_hat + 0.1

  rss <- stance:::compute_residual_sum_squares(Y, S_gamma, U, V, H_v, hrf_basis)
  expect_equal(rss, sum((Y - Y_hat)^2))

  sigma2 <- stance:::update_noise_variance(Y, S_gamma, U, V, H_v, hrf_basis)
  expect_equal(sigma2, rss / (nrow(Y) * ncol(Y)))

  sigma2_prior <- stance:::update_noise_variance(
    Y, S_gamma, U, V, H_v, hrf_basis,
    prior_sigma2 = list(shape = 3, scale = 0.5)
  )
  expect_equal(sigma2_prior,
               (rss + 0.5) / (nrow(Y) * ncol(Y) + 3 + 1))
})
