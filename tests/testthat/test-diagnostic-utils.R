library(testthat)
library(stance)

# Test for convergence diagnostics utility

test_that("check_convergence_diagnostics computes delta correctly", {
  sim <- simulate_fmri_data(V = 6, T = 12, K = 2, algorithm = "CBD", verbose = FALSE)
  cbd <- ContinuousBayesianDecoder$new(Y = sim$Y, K = 2, r = 2)
  cbd$fit(max_iter = 3, verbose = FALSE)
  diag <- check_convergence_diagnostics(cbd)
  elbo <- cbd$get_convergence()$elbo_history
  expected <- abs(elbo[length(elbo)] - elbo[length(elbo) - 1]) / abs(elbo[length(elbo) - 1])
  expect_named(diag, c("converged", "iterations", "delta"))
  expect_equal(diag$delta, expected)
})

# Test for parameter selection helper

test_that("select_best_parameters picks highest scoring combination", {
  cv_results <- data.frame(
    rank = c(2,2,3,3,2,3),
    lambda_h = c(0.1,0.1,0.1,0.1,1,1),
    score = c(0.4,0.6,0.8,1.0,0.7,1.2)
  )
  best <- stance:::select_best_parameters(cv_results)
  expect_equal(best$rank, 3)
  expect_equal(best$lambda_h, 1)
})

# Test for recommendation text helper

test_that("recommend_parameters returns informative message", {
  Y <- matrix(0, nrow = 5, ncol = 8)
  best <- data.frame(rank = 3, lambda_h = 0.5)
  rec <- stance:::recommend_parameters(Y, best)
  expect_match(rec, "Recommended rank 3 and lambda_h 0.5")
  expect_match(rec, "dataset of size 5 x 8")
})

test_that("generate_diagnostic_report writes html output", {
  diag_list <- list(a = 1, b = 2)
  out_dir <- tempfile("diag")
  path <- stance:::generate_diagnostic_report(diag_list, out_dir)
  expect_true(file.exists(path))
  content <- readLines(path)
  expect_true(any(grepl("<h1>Model Diagnostics</h1>", content)))
  expect_true(any(grepl("a", content)))
})

