library(testthat)
library(stance)

# S3-T14 Comprehensive Integration Testing

test_that("test_sprint3_integration runs and returns results", {
  res <- test_sprint3_integration()
  expect_type(res, "list")
  expect_true(all(c("sprint1_baseline", "sprint2_rcpp", "sprint3_gmrf", "sprint3_full") %in% names(res)))
  for (r in res) {
    expect_true(is.numeric(r$time))
    expect_true(is.logical(r$converged))
  }
})
