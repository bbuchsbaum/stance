library(stance)

test_that("prox_tv_dual handles length-2 input", {
  skip_if_not(exists("prox_tv_dual"), "Rcpp functions not compiled")
  x <- c(1, 2)
  lambda <- 0.5
  result <- stance:::prox_tv_dual(x, lambda)
  # prox_tv_dual returns a matrix
  expect_equal(as.vector(result), x)
})
