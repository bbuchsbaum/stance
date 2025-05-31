library(stance)

test_that("simulate_fmri_data validates integer inputs", {
  expect_error(simulate_fmri_data(V = -10, T = 50, K = 2), "V must be a positive integer")
  expect_error(simulate_fmri_data(V = 10, T = 0, K = 2), "T must be a positive integer")
  expect_error(simulate_fmri_data(V = 10, T = 50, K = 0), "K must be a positive integer")
  expect_error(simulate_fmri_data(V = 10, T = 50, K = 2, rank = -1), "rank must be a positive integer")
})

test_that("simulate_fmri_data validates dims argument", {
  expect_error(
    simulate_fmri_data(V = 27, T = 20, K = 2, return_neuroim = TRUE, dims = c(3, 3)),
    "dims must have length 3"
  )
  expect_error(
    simulate_fmri_data(V = 27, T = 20, K = 2, return_neuroim = TRUE, dims = c(3, 3, 4)),
    "product of dims must equal V"
  )
})
