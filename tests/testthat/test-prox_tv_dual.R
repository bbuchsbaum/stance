skip_if_not(exists("prox_tv_dual"), "Rcpp functions not compiled")

# Input validation tests for prox_tv_dual

test_that("prox_tv_dual validates inputs", {
  # Non-finite values in x
  x_bad <- c(1, NA, 3)
  expect_error(
    stance:::prox_tv_dual(x_bad, lambda = 0.5),
    "NaN or Inf" 
  )

  # Negative lambda
  expect_error(
    stance:::prox_tv_dual(c(1, 2, 3), lambda = -0.1),
    "lambda must be non-negative and finite"
  )

  # Non-finite lambda
  expect_error(
    stance:::prox_tv_dual(c(1, 2, 3), lambda = Inf),
    "lambda must be non-negative and finite"
  )

  # Invalid max_iter
  expect_error(
    stance:::prox_tv_dual(c(1, 2, 3), lambda = 0.1, max_iter = 0),
    "max_iter must be positive"
  )

  # Invalid tolerance
  expect_error(
    stance:::prox_tv_dual(c(1, 2, 3), lambda = 0.1, tol = 0),
    "tol must be positive and finite"
  )
})
