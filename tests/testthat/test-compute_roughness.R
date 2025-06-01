library(testthat)
library(stance)

set.seed(123)
H <- matrix(rnorm(20), nrow = 5, ncol = 4)
L <- diag(5) - matrix(1/5, 5, 5)

compute_ref <- function(H, L) {
  total <- 0
  for (b in seq_len(ncol(H))) {
    h <- H[, b]
    total <- total + as.numeric(t(h) %*% L %*% h)
  }
  total / ncol(H)
}

ref_val <- compute_ref(H, L)
fast_val <- compute_roughness(H, L)

test_that("compute_roughness matches loop implementation", {
  expect_equal(fast_val, ref_val)
})
