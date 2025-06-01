library(testthat)
library(stance)


test_that("optimize_convolution_strategy returns sensible defaults", {
  res <- optimize_convolution_strategy(200, 20)
  expect_true(is.list(res))
  expect_true(res$fft >= res$batch_fft)
})


