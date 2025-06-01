library(testthat)
library(stance)

skip_if_not_installed("neuroim2")

# Integration test ensuring GMRF prior enforces spatial smoothness

test_that("GMRF prior smooths voxel HRFs", {
  sim <- simulate_fmri_data(
    V = 1000, T = 100, K = 3,
    algorithm = "CBD",
    spatial_hrf = list(correlation_length = 5, smoothness = 0.8),
    dims = c(10,10,10),
    verbose = FALSE
  )

  mask_data <- array(TRUE, c(10,10,10))
  mask <- neuroim2::NeuroVol(mask_data,
                             neuroim2::NeuroSpace(dim = c(10,10,10)))

  cbd_no <- ContinuousBayesianDecoder$new(
    Y = sim$Y, K = 3, r = 10,
    use_gmrf = FALSE,
    hrf_basis = "canonical"
  )

  cbd_gm <- ContinuousBayesianDecoder$new(
    Y = sim$Y, K = 3, r = 10,
    use_gmrf = TRUE, lambda_h = 10,
    mask = mask,
    hrf_basis = "canonical"
  )

  cbd_no$fit(max_iter = 5, verbose = FALSE)
  cbd_gm$fit(max_iter = 5, verbose = FALSE)

  L <- create_gmrf_laplacian_neuroim2(neuroim2::NeuroSpace(dim = c(10,10,10)))
  rough_no <- compute_roughness(cbd_no$get_hrf_estimates(), L)
  rough_gm <- compute_roughness(cbd_gm$get_hrf_estimates(), L)

  expect_true(rough_gm < 0.7 * rough_no)
})
