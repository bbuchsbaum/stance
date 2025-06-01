library(testthat)
library(stance)

# Use R implementations to avoid compiled code during tests
formals(stance:::ContinuousBayesianDecoder$public_methods$initialize)$engine <- "R"

test_check("stance")
