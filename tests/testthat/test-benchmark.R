library(testthat)
library(stance)

# S3-T16 Performance Benchmarking

test_that("benchmark_against_proposal runs and returns expected structure", {
  res <- benchmark_against_proposal(max_iter = 1)
  expect_type(res, "list")
  expect_true(all(c("roi", "parcellated", "full_wb") %in% names(res)))
  for (r in res) {
    expect_true(is.numeric(r$time_minutes))
    expect_true(is.numeric(r$memory_gb))
    expect_true(is.logical(r$meets_target))
    expect_true(is.numeric(r$efficiency))
  }
})
