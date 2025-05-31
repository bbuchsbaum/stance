library(testthat)
library(stance)

test_that("create_chain_laplacian returns valid Laplacian", {
  L <- create_chain_laplacian(5)
  expect_s4_class(L, "dgCMatrix")
  expect_equal(dim(L), c(5, 5))
  expect_true(Matrix::isSymmetric(L))
  expect_equal(Matrix::rowSums(L), rep(0, 5))
  expect_equal(diag(L), c(1, rep(2, 3), 1))
})

test_that("create_gmrf_laplacian_neuroim2 produces correct structure", {
  skip_if_not_installed("neuroim2")
  space <- neuroim2::NeuroSpace(dim = c(2, 2, 1), spacing = c(1, 1, 1))
  L <- create_gmrf_laplacian_neuroim2(space)
  expect_s4_class(L, "dgCMatrix")
  expect_equal(dim(L), c(4, 4))
  expect_true(Matrix::isSymmetric(L))
  expect_equal(Matrix::rowSums(L), rep(0, 4))
  expect_equal(diag(L), rep(2, 4))
})
