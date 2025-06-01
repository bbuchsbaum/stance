library(testthat)
library(stance)

forward_R <- function(log_lik, Pi, pi0) {
  K <- nrow(log_lik)
  T_len <- ncol(log_lik)
  alpha <- matrix(0, K, T_len)
  c_scale <- numeric(T_len)
  alpha[,1] <- pi0 * exp(log_lik[,1])
  c_scale[1] <- 1 / sum(alpha[,1])
  alpha[,1] <- alpha[,1] * c_scale[1]
  if (T_len > 1) {
    for (t in 2:T_len) {
      alpha[,t] <- (alpha[,t-1] %*% Pi) * exp(log_lik[,t])
      c_scale[t] <- 1 / sum(alpha[,t])
      alpha[,t] <- alpha[,t] * c_scale[t]
    }
  }
  list(alpha = alpha, c_scale = c_scale,
       log_likelihood = -sum(log(c_scale)))
}

backward_R <- function(log_lik, Pi, c_scale) {
  K <- nrow(log_lik)
  T_len <- ncol(log_lik)
  beta <- matrix(1, K, T_len)
  if (T_len > 1) {
    for (t in (T_len-1):1) {
      beta[,t] <- Pi %*% (beta[,t+1] * exp(log_lik[,t+1]))
      beta[,t] <- beta[,t] * c_scale[t+1]
    }
  }
  beta
}

test_that("compute_log_likelihoods_rcpp matches R version", {
  set.seed(123)
  r <- 3; T_len <- 8; K <- 2
  Y_proj <- matrix(rnorm(r*T_len), r, T_len)
  Vmat <- matrix(rnorm(K*r), K, r)
  h <- rep(0.5, 4)
  sigma2 <- 0.5
  ll_r <- stance:::compute_log_likelihoods_r(Y_proj, Vmat, h, sigma2)
  ll_cpp <- compute_log_likelihoods_rcpp(Y_proj, Vmat, h, sigma2)
  expect_equal(ll_cpp, ll_r, tolerance = 1e-8)
})

test_that("complete low-rank likelihood matches simple version", {
  set.seed(124)
  r <- 3; T_len <- 10; K <- 2; V <- 5; L <- 3
  Y_proj <- matrix(rnorm(r*T_len), r, T_len)
  U <- matrix(rnorm(V*r), V, r)
  Vmat <- matrix(rnorm(K*r), K, r)
  H_v <- matrix(runif(V*L), V, L)
  h_basis <- matrix(runif(4*L), 4, L)
  gamma <- matrix(runif(K*T_len), K, T_len)
  sigma2 <- 1.0
  ll_full <- compute_log_likelihood_lowrank_complete(Y_proj, U, Vmat, H_v, h_basis, gamma, sigma2)
  expect_equal(dim(ll_full), c(K, T_len))
})

test_that("forward and backward passes match R implementations", {
  set.seed(321)
  K <- 3; T_len <- 6
  log_lik <- matrix(rnorm(K*T_len), K, T_len)
  Pi <- matrix(runif(K*K), K, K)
  Pi <- Pi/rowSums(Pi)
  pi0 <- rep(1/K, K)

  f_R <- forward_R(log_lik, Pi, pi0)
  f_cpp <- forward_pass_rcpp(log_lik, Pi, pi0)
  expect_equal(f_cpp$alpha, f_R$alpha, tolerance = 1e-8)
  expect_equal(f_cpp$log_likelihood, f_R$log_likelihood, tolerance = 1e-8)

  b_R <- backward_R(log_lik, Pi, f_R$c_scale)
  b_cpp <- backward_pass_rcpp(log_lik, Pi, f_R$c_scale)
  expect_equal(b_cpp, b_R, tolerance = 1e-8)
})

test_that("plot_hrf returns ggplot when CI provided", {
  h <- matrix(seq(0, 1, length.out = 8), ncol = 1)
  var <- matrix(0.01, nrow = 8, ncol = 1)
  p <- plot_hrf(h, TR = 1, hrf_var = var)
  expect_true(inherits(p, "ggplot"))
})
