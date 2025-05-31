# Verify compute_tv_rcpp output

test_that("compute_tv_rcpp matches manual calculation", {
  skip_if_not(exists("compute_tv_rcpp"), "Rcpp functions not compiled")

  X <- matrix(c(1, 3, 2, 5, 4, 6), nrow = 2, byrow = TRUE)
  tv_manual <- sum(abs(X[, -1] - X[, -ncol(X)]))
  tv_cpp <- stance:::compute_tv_rcpp(X)

  expect_equal(tv_cpp, tv_manual)
})
