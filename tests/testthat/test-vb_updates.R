library(testthat)
library(stance)

# Synthetic data for VB update tests
set.seed(123)
V <- 5
T_obs <- 4
K <- 2
Y <- matrix(rnorm(V * T_obs), V, T_obs)
U <- matrix(rnorm(V * 2), V, 2)
V_mat <- matrix(rnorm(K * 2), K, 2)
H_v <- matrix(rnorm(V), V, 1)
hrf_basis <- matrix(1, 1, 1)
Pi <- matrix(c(0.7, 0.3,
               0.4, 0.6),
             nrow = K, byrow = TRUE)
pi0 <- c(0.6, 0.4)
params <- list(
  U = U,
  V = V_mat,
  H_v = H_v,
  hrf_basis = hrf_basis,
  Pi = Pi,
  pi0 = pi0,
  sigma2 = 1,
  L_gmrf = diag(V),
  lambda_H_prior = 1,
  sigma2_prior = NULL,
  prior_Pi = NULL,
  prior_pi0 = NULL
)
config <- list(engine = "cpp")

# Stub implementations to avoid missing compiled code
fb_stub <- function(Y, U, V, H_v, hrf_basis, Pi, pi0, sigma2, engine) {
  K <- nrow(Pi)
  Tt <- ncol(Y)
  list(
    gamma = matrix(1 / K, K, Tt),
    xi = array(1 / (K * K), c(K, K, Tt - 1)),
    log_likelihood = -1
  )
}
spatial_stub <- function(Y, S_gamma, H_v, hrf_basis, current_U, current_V) {
  list(U = current_U, V = current_V)
}
hrf_stub <- function(Y, S_gamma, U, V, hrf_basis, L_gmrf, lambda_H_prior, sigma2) {
  H_v
}

assignInNamespace("forward_backward_algorithm", fb_stub, ns = "stance")
assignInNamespace("update_spatial_components_cpp", spatial_stub, ns = "stance")
assignInNamespace("update_hrf_coefficients_gmrf_cpp", hrf_stub, ns = "stance")

test_that("vb_e_step returns valid dimensions and probabilities", {
  res <- stance:::vb_e_step(Y, params, config)
  expect_equal(dim(res$S_gamma), c(K, T_obs))
  expect_equal(dim(res$S_xi), c(K, K, T_obs - 1))
  expect_equal(apply(res$S_gamma, 2, sum), rep(1, T_obs))
  expect_equal(apply(res$S_xi, 3, sum), rep(1, T_obs - 1))
})

vb_params <- stance:::vb_e_step(Y, params, config)

test_that("vb_m_step updates return correct structure", {
  out <- stance:::vb_m_step(Y, vb_params, params, config)
  expect_equal(dim(out$U), dim(U))
  expect_equal(dim(out$V), dim(V_mat))
  expect_equal(dim(out$H_v), dim(H_v))
  expect_equal(dim(out$Pi), dim(Pi))
  expect_equal(length(out$pi0), length(pi0))
  expect_true(is.numeric(out$sigma2))
})

test_that("compute_elbo produces a numeric scalar", {
  elbo <- stance:::compute_elbo(Y, vb_params, params, config)
  expect_type(elbo, "double")
  expect_length(elbo, 1)
})

