test_that("validate_fmri_input works with matrix input", {
  # Create test matrix
  Y_mat <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
  
  # Test basic validation
  result <- validate_fmri_input(Y_mat)
  expect_equal(result$type, "matrix")
  expect_equal(result$dims["V"], 100)
  expect_equal(result$dims["T"], 50)
  expect_true(is.matrix(result$data))
  expect_null(result$space)
  expect_null(result$mask)
  
  # Test with expected dimensions
  result2 <- validate_fmri_input(Y_mat, expected_dims = list(V = 100, T = 50))
  expect_equal(result2$dims["V"], 100)
  
  # Test dimension mismatch
  expect_error(
    validate_fmri_input(Y_mat, expected_dims = list(V = 200)),
    "Expected 200 voxels but got 100"
  )
  
  # Test non-finite values
  Y_mat_bad <- Y_mat
  Y_mat_bad[1, 1] <- NA
  expect_error(
    validate_fmri_input(Y_mat_bad),
    "Data contains non-finite values"
  )
})

test_that("validate_fmri_input handles data.frame input", {
  # Create test data.frame
  Y_df <- as.data.frame(matrix(rnorm(100 * 50), nrow = 100, ncol = 50))

  result <- validate_fmri_input(Y_df)
  expect_equal(result$type, "data.frame")
  expect_true(is.matrix(result$data))
  expect_true(is.numeric(result$data))
  expect_equal(dim(result$data), c(100, 50))
})

test_that("validate_fmri_input errors with non-numeric data.frame", {
  Y_bad <- data.frame(a = rnorm(10), b = letters[1:10])
  expect_error(
    validate_fmri_input(Y_bad),
    "non-numeric columns"
  )
})

test_that("extract_data_matrix handles different input types", {
  # Test with matrix
  mat <- matrix(1:20, nrow = 4, ncol = 5)
  result <- extract_data_matrix(mat)
  expect_equal(result$data, mat)
  expect_equal(result$metadata$class, "matrix")
  
  # Test preserve_attributes = FALSE
  result2 <- extract_data_matrix(mat, preserve_attributes = FALSE)
  expect_equal(length(result2$metadata), 0)
})

test_that("check_temporal_alignment validates design matrices", {
  # Create test data
  Y_info <- list(dims = c(V = 100, T = 200), space = NULL)
  design_mat <- matrix(0, nrow = 3, ncol = 200)
  TR <- 2.0
  
  # Should pass
  expect_true(check_temporal_alignment(Y_info, design_mat, TR))
  
  # Wrong dimensions
  bad_design <- matrix(0, nrow = 3, ncol = 150)
  expect_error(
    check_temporal_alignment(Y_info, bad_design, TR),
    "Design matrix has 150 columns but data has 200 timepoints"
  )
  
  # Test with event data.frame
  events <- data.frame(
    onset = c(0, 10, 20),
    duration = c(2, 2, 2)
  )
  expect_true(check_temporal_alignment(Y_info, events, TR))
  
  # Events extend beyond data
  bad_events <- data.frame(
    onset = c(0, 10, 500),
    duration = c(2, 2, 2)
  )
  expect_error(
    check_temporal_alignment(Y_info, bad_events, TR),
    "Events extend to .* but data only covers"
  )
})

test_that("create_output_structure works for both algorithms", {
  # Test data
  spatial_maps <- matrix(rnorm(100 * 3), nrow = 100, ncol = 3)
  temporal_activations <- matrix(rnorm(3 * 50), nrow = 3, ncol = 50)
  metadata <- list(source = "test")
  
  # Test CLD output
  cld_out <- create_output_structure(
    spatial_maps, temporal_activations, metadata, "CLD",
    additional_info = list(
      objective_values = c(100, 90, 80),
      convergence_info = list(iterations = 50)
    )
  )
  
  expect_s3_class(cld_out, "CLD_output")
  expect_s3_class(cld_out, "stance_output")
  expect_equal(cld_out$algorithm, "CLD")
  expect_equal(cld_out$dims$V, 100)
  expect_equal(cld_out$dims$K, 3)
  expect_equal(cld_out$dims$T, 50)
  expect_true(!is.null(cld_out$objective_values))
  
  # Test CBD output
  cbd_out <- create_output_structure(
    spatial_maps, temporal_activations, metadata, "CBD",
    additional_info = list(
      uncertainty = matrix(0.1, nrow = 3, ncol = 50),
      posterior_params = list(alpha = 1, beta = 2)
    )
  )
  
  expect_s3_class(cbd_out, "CBD_output")
  expect_equal(cbd_out$algorithm, "CBD")
  expect_true(!is.null(cbd_out$uncertainty))
  expect_true(!is.null(cbd_out$posterior_params))
})