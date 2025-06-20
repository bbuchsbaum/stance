library(testthat)
library(stance)

# Basic initialization test for ContinuousBayesianDecoder

test_that("ContinuousBayesianDecoder initializes on simulated data", {
  # Simulate small dataset
  V <- 20
  Tlen <- 30
  K <- 3
  sim <- simulate_fmri_data(V = V, T = Tlen, K = K, algorithm = "CBD", verbose = FALSE)

  # Initialize decoder
  cbd <- ContinuousBayesianDecoder$new(
    Y = sim$Y,
    K = K,
    r = 5,
    hrf_basis = "canonical",
    engine = "R"
  )

  expect_s3_class(cbd, "ContinuousBayesianDecoder")

  priv <- cbd$.__enclos_env__$private

  expect_equal(priv$.n_voxels, V)
  expect_equal(priv$.T, Tlen)
  expect_equal(priv$.K, K)
  expect_equal(priv$.r, 5)

  expect_equal(dim(priv$.U), c(V, 5))
  expect_equal(dim(priv$.V), c(K, 5))
  expect_equal(dim(priv$.S_gamma), c(K, Tlen))
  expect_equal(dim(priv$.S_xi), c(K, K, Tlen - 1))
  expect_true(is.numeric(priv$.sigma2))

  expect_true(is.matrix(priv$.hrf_basis))
  expect_equal(dim(priv$.H_v), c(V, ncol(priv$.hrf_basis)))
  expect_gt(sum(abs(priv$.H_v)), 0)
})
