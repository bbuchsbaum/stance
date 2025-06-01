library(testthat)
library(stance)

context("GMRF batched solver")

test_that("solve_gmrf_batched returns correct dimensions", {
  V <- 10
  L <- 2
  XtX <- Matrix::Diagonal(V)
  Lmat <- create_chain_laplacian(V)
  XtY <- matrix(rnorm(V * L), V, L)
  H <- stance:::solve_gmrf_batched(XtX, Lmat, XtY, lambda_h = 0.5, block_size = 4L)
  expect_equal(dim(H), c(V, L))

  Q <- XtX + 0.5 * Lmat
  H_ref <- Matrix::solve(Q, XtY)
  expect_equal(as.numeric(H[1,1]), as.numeric(H_ref[1,1]), tolerance = 1e-6)
})
