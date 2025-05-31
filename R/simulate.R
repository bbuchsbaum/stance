#' Simulation Utilities
#'
#' Generate synthetic fMRI data for testing.
#'
#' @name simulate
NULL

#' Simulate fMRI Data
#'
#' @param V Number of voxels
#' @param T Number of time points
#' @param K Number of states
#' @return List with Y matrix and design
#' @export
simulate_fmri_data <- function(V = 100, T = 200, K = 3) {
  W <- matrix(rnorm(V * K), V, K)
  S <- matrix(rbinom(K * T, 1, 0.5), K, T)
  h <- setup_hrf_kernel()
  X <- convolve_with_hrf(S, h)
  Y <- W %*% X + matrix(rnorm(V * T, sd = 0.1), V, T)
  list(Y = Y, S = S)
}
