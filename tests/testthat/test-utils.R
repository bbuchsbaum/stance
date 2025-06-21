library(stance)

# tests for null coalescing operator

test_that("%||% returns first value when not NULL", {
  expect_equal(1 %||% 2, 1)
  expect_equal("a" %||% "b", "a")
})

test_that("%||% returns second value when first is NULL", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(NULL %||% NULL %||% 7, 7)
})

# tests for normalize_rows utility

test_that("normalize_rows scales rows to sum to one", {
  X <- matrix(c(1,2,3,4), nrow = 2, byrow = TRUE)
  res <- stance:::normalize_rows(X)
  expect_equal(rowSums(res), c(1,1))
  expect_equal(res[1,], c(1/3,2/3), tolerance = 1e-8)
})

# zero rows should stay zero without NaNs

test_that("normalize_rows handles zero rows", {
  X <- matrix(c(0,0,1,1), nrow = 2, byrow = TRUE)
  res <- stance:::normalize_rows(X)
  expect_false(any(is.nan(res)))
  expect_equal(rowSums(res), c(0,1))
})

# tests for time parsing helpers

test_that("parse_time parses ranges with minutes and hours", {
  expect_equal(stance:::parse_time("5-15 minutes"), 10)
  expect_equal(stance:::parse_time("4-12 hours"), 480)
})

test_that("parse_minutes is a thin wrapper for parse_time", {
  expect_equal(stance:::parse_minutes("90 minutes"), 90)
  expect_equal(stance:::parse_minutes("2-4 hours"), 180)
})

