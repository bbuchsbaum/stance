library(microbenchmark)
library(stance)

# Naive loop implementations for comparison
compute_data_likelihood_term_loop <- function(Y, S_gamma, U, V, sigma2) {
  V_voxels <- nrow(Y)
  T <- ncol(Y)
  K <- nrow(S_gamma)
  log_lik <- -0.5 * V_voxels * T * log(2 * pi * sigma2)
  reconstruction_error <- 0
  for (t in seq_len(T)) {
    y_t <- Y[, t]
    y_hat_t <- numeric(V_voxels)
    for (k in seq_len(K)) {
      gamma_kt <- S_gamma[k, t]
      if (gamma_kt > 1e-10) {
        w_k_contribution <- U %*% V[k, ]
        y_hat_t <- y_hat_t + gamma_kt * w_k_contribution
      }
    }
    reconstruction_error <- reconstruction_error + sum((y_t - y_hat_t)^2)
  }
  log_lik - 0.5 * reconstruction_error / sigma2
}

compute_residual_sum_squares_loop <- function(Y, S_gamma, U, V) {
  V_voxels <- nrow(Y)
  T <- ncol(Y)
  K <- nrow(S_gamma)
  rss <- 0
  for (t in seq_len(T)) {
    y_t <- Y[, t]
    y_hat_t <- numeric(V_voxels)
    for (k in seq_len(K)) {
      gamma_kt <- S_gamma[k, t]
      if (gamma_kt > 1e-10) {
        w_k <- U %*% V[k, ]
        y_hat_t <- y_hat_t + gamma_kt * w_k
      }
    }
    rss <- rss + sum((y_t - y_hat_t)^2)
  }
  rss
}

# Set up small synthetic data
set.seed(123)
V <- 200
Tt <- 50
K <- 4
rank <- 3
sim <- simulate_fmri_data(V, Tt, K, rank, algorithm = "CBD", verbose = FALSE)
U <- sim$U
Vmat <- sim$V_mat
Y <- sim$Y
Sg <- sim$S
sigma2 <- 1

bench <- microbenchmark(
  vector_lik = compute_data_likelihood_term(Y, Sg, U, Vmat, NULL, NULL, sigma2),
  loop_lik    = compute_data_likelihood_term_loop(Y, Sg, U, Vmat, sigma2),
  vector_rss  = compute_residual_sum_squares(Y, Sg, U, Vmat, NULL, NULL),
  loop_rss    = compute_residual_sum_squares_loop(Y, Sg, U, Vmat),
  times = 10
)
print(bench)
