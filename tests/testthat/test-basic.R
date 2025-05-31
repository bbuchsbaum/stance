library(stance)

test_that("simulation works", {
  sim <- simulate_fmri_data(10, 20, 2)
  expect_equal(dim(sim$Y), c(10, 20))
  expect_equal(dim(sim$S), c(2, 20))
})
