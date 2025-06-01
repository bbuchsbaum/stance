library(stance)

test_that("fista_tv_rcpp multi-threading matches single-thread", {
  skip_if_not(exists("fista_tv_rcpp"), "Rcpp functions not compiled")
  openmp_info <- check_openmp_support()
  if (openmp_info$available) {
    set.seed(123)
    V <- 20; T_len <- 30; K <- 3
    Y <- matrix(rnorm(V * T_len), V, T_len)
    W <- matrix(rnorm(V * K), V, K)
    hrf <- rep(1/3, 3)
    WtY <- t(W) %*% Y
    L <- estimate_lipschitz_rcpp(W, hrf)
    X_init <- matrix(0, K, T_len)
    res1 <- fista_tv_rcpp(WtY, W, hrf, 0.01, L, X_init,
                          max_iter = 5L, tol = 1e-4,
                          verbose = FALSE, n_threads = 1L)
    res2 <- fista_tv_rcpp(WtY, W, hrf, 0.01, L, X_init,
                          max_iter = 5L, tol = 1e-4,
                          verbose = FALSE, n_threads = 2L)
    expect_equal(res1$X_hat, res2$X_hat, tolerance = 1e-6)
  }
})

test_that("compute_gradient_fista_precomp_rcpp multi-threading", {
  skip_if_not(exists("compute_gradient_fista_precomp_rcpp"), "Rcpp functions not compiled")
  openmp_info <- check_openmp_support()
  if (openmp_info$available) {
    set.seed(456)
    K <- 2; T_len <- 20; V <- 5
    W <- matrix(rnorm(V * K), V, K)
    X <- matrix(rnorm(K * T_len), K, T_len)
    hrf <- c(0.2, 0.5, 0.3)
    H_star_X <- convolve_rows_rcpp(X, hrf, n_threads = 1L)
    WtW <- t(W) %*% W
    WtY <- WtW %*% H_star_X
    grad1 <- compute_gradient_fista_precomp_rcpp(WtY, WtW, H_star_X, hrf, n_threads = 1L)
    grad2 <- compute_gradient_fista_precomp_rcpp(WtY, WtW, H_star_X, hrf, n_threads = 2L)
    expect_equal(grad1, grad2, tolerance = 1e-8)
  }
})
