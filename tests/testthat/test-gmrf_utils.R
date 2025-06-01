library(testthat)
library(stance)

test_that("create_chain_laplacian returns valid Laplacian", {
  L <- create_chain_laplacian(5)
  expect_s4_class(L, "dgCMatrix")
  expect_equal(dim(L), c(5, 5))
  expect_true(Matrix::isSymmetric(L))
  expect_equal(Matrix::rowSums(L), rep(0, 5))
  expect_equal(as.vector(Matrix::diag(L)), c(1, rep(2, 3), 1))
})

test_that("create_gmrf_laplacian_neuroim2 produces correct structure", {
  skip_if_not_installed("neuroim2")
  space <- neuroim2::NeuroSpace(dim = c(2, 2, 1), spacing = c(1, 1, 1))
  L <- create_gmrf_laplacian_neuroim2(space)
  expect_s4_class(L, "dgCMatrix")
  expect_equal(dim(L), c(4, 4))
  expect_true(Matrix::isSymmetric(L))
  expect_equal(Matrix::rowSums(L), rep(0, 4))
  expect_equal(as.vector(Matrix::diag(L)), rep(2, 4))
})

test_that("get_spatial_neighbors works with NeuroVol mask", {
  skip_if_not_installed("neuroim2")
  mask_arr <- array(1, dim = c(2,2,1))
  vol <- neuroim2::NeuroVol(mask_arr, neuroim2::NeuroSpace(c(2,2,1)))
  info <- get_spatial_neighbors(vol)
  expect_equal(length(info$neighbors), 4)
  expect_equal(info$neighbors,
               list(c(2L,3L), c(1L,4L), c(4L,1L), c(3L,2L)))
  L <- create_gmrf_laplacian(info$neighbors, length(info$neighbors))
  expect_true(Matrix::isSymmetric(L))
  expect_equal(Matrix::rowSums(L), rep(0, 4))
})

test_that("get_spatial_neighbors handles logical mask", {
  mask <- array(TRUE, dim = c(2,1,1))
  info <- get_spatial_neighbors(mask)
  L <- create_gmrf_laplacian(info$neighbors, length(info$neighbors))
  expect_equal(dim(L), c(2,2))
  expect_equal(info$neighbors, list(2L,1L))
  expect_true(Matrix::isSymmetric(L))
  expect_equal(Matrix::rowSums(L), rep(0, 2))
})
