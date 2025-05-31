test_that("vectorized spatial smoothing matches legacy implementation", {
  arr <- array(rnorm(1000), dim = c(10, 10, 10))
  old <- stance:::spatial_smooth_3d_loop(arr)
  new <- stance:::spatial_smooth_3d(arr)
  expect_equal(new, old, tolerance = 1e-8)
})
