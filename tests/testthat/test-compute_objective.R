library(testthat)
library(stance)

test_that("compute_objective matches manual calculation", {
  W <- diag(2)
  X <- matrix(c(1, 2, 3,
                4, 5, 6), nrow = 2, byrow = TRUE)
  hrf <- c(0.5, 0.5)
  HX <- convolve_with_hrf(X, hrf)
  Y <- W %*% HX
  lambda <- 0.3

  manual <- 0.5 * sum((Y - W %*% HX)^2) +
    lambda * stance:::compute_tv(X)
  expect_equal(stance:::compute_objective(Y, W, X, hrf, lambda), manual)

  Y2 <- Y + 0.1
  manual2 <- 0.5 * sum((Y2 - W %*% HX)^2) +
    lambda * stance:::compute_tv(X)
  val2 <- stance:::compute_objective(Y2, W, X, hrf, lambda)
  expect_equal(val2, manual2)
  expect_gt(val2, manual)
})
