library(stance)

# fast_matmul caching behavior

test_that("fast_matmul caches results", {
  skip("fast_matmul not yet implemented")
  A <- matrix(1:4, 2, 2)
  B <- diag(2)
  env <- new.env()

  result1 <- stance:::fast_matmul(A, B, cache_key = "prod", cache_env = env)
  # change B but use same cache key
  B2 <- B * 2
  result2 <- stance:::fast_matmul(A, B2, cache_key = "prod", cache_env = env)

  # should return cached result, not recomputed with B2
  expect_equal(result2, result1)
  # cached value stored in environment
  expect_equal(env$prod, result1)
})

# process_blocks execution and assembly

test_that("process_blocks applies function by block", {
  skip("process_blocks not yet implemented")
  X <- matrix(1:20, nrow = 5, ncol = 4)
  counter <- new.env(); counter$n <- 0
  fun <- function(block) {
    counter$n <- counter$n + 1
    block * 2
  }

  res <- stance:::process_blocks(X, fun, block_size = 2)
  expect_equal(res, X * 2)
  # there should be three blocks: 2 + 2 + 1 rows
  expect_equal(counter$n, 3)
})

# estimate_operator_norm convergence

test_that("estimate_operator_norm converges on simple matrix", {
  A <- matrix(c(2, 0, 0, 1), 2, 2)
  est <- stance:::estimate_operator_norm(A, dim = 2, max_iter = 100, tol = 1e-8)
  expect_equal(est, norm(A, type = "2"), tolerance = 1e-4)
})

# project_orthogonal orthonormal columns

test_that("project_orthogonal returns orthonormal columns", {
  set.seed(1)
  U <- matrix(rnorm(20), 5, 4)
  Uo <- stance:::project_orthogonal(U)
  prod <- crossprod(Uo)
  expect_equal(dim(Uo), c(5, 4))
  expect_equal(prod, diag(4), tolerance = 1e-10)
})

# create_progress_bar output

test_that("create_progress_bar produces expected text", {
  pb <- stance:::create_progress_bar(total = 3, width = 10)
  out <- capture.output({
    pb(1); pb(2); pb(3)
  })
  expect_true(any(grepl("Progress: |", out, fixed = TRUE)))
  expect_true(any(grepl("100%", out)))
})

