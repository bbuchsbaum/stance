library(testthat)
library(stance)

set.seed(2)
sim <- simulate_fmri_data(V = 8, T = 16, K = 2, algorithm = "CBD", verbose = FALSE)

res <- cbd_cross_validate(sim$Y, rank_values = c(2, 3),
                          lambda_h_values = c(0.1, 1), n_folds = 2,
                          use_parallel = FALSE)

test_that("cbd_cross_validate returns expected structure", {
  expect_true(is.list(res))
  expect_true(all(c("results", "best_params", "recommendation") %in% names(res)))
  expect_equal(nrow(res$results), 2*2*2)
})
