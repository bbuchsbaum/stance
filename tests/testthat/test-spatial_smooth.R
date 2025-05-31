test_that("spatial smoothing using neuroim2 works", {
  arr <- array(rnorm(1000), dim = c(10, 10, 10))
  
  # Test that the function runs without error
  smoothed <- stance:::spatial_smooth_3d(arr, sigma = 1)
  
  # Check output dimensions match input
  expect_equal(dim(smoothed), dim(arr))
  
  # Check that smoothing reduced variance (as expected)
  expect_true(var(as.vector(smoothed)) < var(as.vector(arr)))
  
  # Check no NAs introduced
  expect_false(any(is.na(smoothed)))
})
