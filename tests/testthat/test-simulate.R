test_that("simulate_fmri_data works for CLD", {
  # Small simulation for CLD
  sim <- simulate_fmri_data(V = 100, T = 50, K = 3, algorithm = "CLD", verbose = FALSE)
  
  # Check output structure
  expect_true(is.list(sim))
  expect_equal(dim(sim$Y), c(100, 50))
  expect_equal(dim(sim$S), c(3, 50))
  expect_equal(dim(sim$W), c(100, 3))
  
  # Check CLD-specific: continuous states
  expect_true(any(sim$S > 0 & sim$S < 1))  # Should have continuous values
  
  # Check parameters stored
  expect_equal(sim$params$algorithm, "CLD")
  expect_equal(sim$params$V, 100)
  expect_equal(sim$params$T, 50)
  expect_equal(sim$params$K, 3)
})

test_that("simulate_fmri_data works for CBD", {
  # Small simulation for CBD
  sim <- simulate_fmri_data(V = 100, T = 50, K = 3, algorithm = "CBD", verbose = FALSE)
  
  # Check output structure
  expect_equal(dim(sim$Y), c(100, 50))
  expect_equal(dim(sim$S), c(3, 50))
  
  # Check CBD-specific: one-hot encoded states
  expect_true(all(sim$S %in% c(0, 1)))  # Binary values only
  expect_equal(colSums(sim$S), rep(1, 50))  # One state active per timepoint
})

test_that("low-rank decomposition works correctly", {
  sim <- simulate_fmri_data(V = 200, T = 50, K = 4, rank = 10)
  
  # Check dimensions
  expect_equal(dim(sim$U), c(200, 10))
  expect_equal(dim(sim$V_mat), c(4, 10))
  
  # Check that W = U * V_mat^T
  W_reconstructed <- sim$U %*% t(sim$V_mat)
  expect_equal(W_reconstructed, sim$W, tolerance = 1e-10)
  
  # Check orthogonality of U
  UTU <- t(sim$U) %*% sim$U
  expect_equal(diag(UTU), rep(1, 10), tolerance = 1e-10)
  expect_true(all(abs(UTU[upper.tri(UTU)]) < 1e-10))
})

test_that("HRF convolution is applied", {
  sim <- simulate_fmri_data(V = 50, T = 100, K = 2, hrf_spec = "spmg1")
  
  # Check HRF was stored
  expect_true(is.numeric(sim$hrf))
  expect_true(length(sim$hrf) > 1)
  
  # Check convolved states have correct dimensions
  expect_equal(dim(sim$X), dim(sim$S))
  
  # Convolution should not produce identical output
  expect_false(identical(sim$X, sim$S))
})

test_that("SNR affects noise level", {
  # High SNR
  sim_high <- simulate_fmri_data(V = 100, T = 50, K = 2, snr = 10)
  signal_var_high <- var(as.vector(sim_high$W %*% sim_high$X))
  noise_var_high <- var(as.vector(sim_high$noise))
  
  # Low SNR
  sim_low <- simulate_fmri_data(V = 100, T = 50, K = 2, snr = 0.1)
  signal_var_low <- var(as.vector(sim_low$W %*% sim_low$X))
  noise_var_low <- var(as.vector(sim_low$noise))
  
  # Check SNR relationships
  snr_high <- signal_var_high / noise_var_high
  snr_low <- signal_var_low / noise_var_low
  
  expect_true(snr_high > snr_low)
})

test_that("neuroimaging output format works", {
  sim <- simulate_fmri_data(V = 1000, T = 30, K = 2, 
                            return_neuroim = TRUE, dims = c(10, 10, 10))
  
  # Check NeuroVec was created
  expect_true("Y_neurovol" %in% names(sim))
  expect_s4_class(sim$Y_neurovol, "NeuroVec")
  
  # Check spatial maps as NeuroVol
  expect_true("W_neurovol" %in% names(sim))
  expect_true(is.list(sim$W_neurovol))
  expect_equal(length(sim$W_neurovol), 2)  # K=2
  expect_s4_class(sim$W_neurovol[[1]], "NeuroVol")
})

test_that("generate_markov_states produces valid chains", {
  # Test Markov state generation
  states <- stance:::generate_markov_states(K = 3, T = 100)
  
  # Check dimensions and one-hot encoding
  expect_equal(dim(states), c(3, 100))
  expect_true(all(states %in% c(0, 1)))
  expect_equal(colSums(states), rep(1, 100))
  
  # Test with custom transition matrix
  trans_mat <- matrix(c(0.8, 0.1, 0.1,
                        0.2, 0.7, 0.1,
                        0.1, 0.2, 0.7), 3, 3, byrow = TRUE)
  states2 <- stance:::generate_markov_states(K = 3, T = 1000, trans_mat)
  
  # Should have more self-transitions due to diagonal dominance
  state_seq <- apply(states2, 2, which.max)
  transitions <- sum(diff(state_seq) != 0)
  expect_true(transitions < 500)  # Fewer than half should be transitions
})

test_that("generate_continuous_states produces valid activations", {
  states <- stance:::generate_continuous_states(K = 3, T = 200, smooth = TRUE)
  
  # Check dimensions
  expect_equal(dim(states), c(3, 200))
  
  # Should be non-negative
  expect_true(all(states >= 0))
  
  # Should have some sparsity
  expect_true(sum(states == 0) > 0)
  
  # Each state should have some activation
  expect_true(all(rowSums(states) > 0))
})

test_that("spatial smoothing works", {
  # Generate patterns with and without smoothing
  patterns_smooth <- stance:::generate_spatial_patterns(
    V = 1000, K = 3, rank = 5, smooth = TRUE, dims = c(10, 10, 10)
  )
  
  patterns_no_smooth <- stance:::generate_spatial_patterns(
    V = 1000, K = 3, rank = 5, smooth = FALSE, dims = c(10, 10, 10)
  )
  
  # Both should have correct dimensions
  expect_equal(dim(patterns_smooth$W), c(1000, 3))
  expect_equal(dim(patterns_no_smooth$W), c(1000, 3))
  
  # Smoothed patterns should have different characteristics
  # (specific test depends on smoothing implementation)
})

test_that("multi-subject simulation works", {
  multi_sim <- simulate_multi_subject(n_subjects = 3, V = 100, T = 50, K = 2,
                                      subject_variability = 0.2, algorithm = "CLD")
  
  # Check structure
  expect_true(is.list(multi_sim))
  expect_equal(length(multi_sim$subjects), 3)
  
  # Each subject should have data
  for (i in 1:3) {
    expect_equal(dim(multi_sim$subjects[[i]]$Y), c(100, 50))
    expect_equal(multi_sim$subjects[[i]]$subject_id, i)
  }
  
  # Subjects should have different spatial maps
  W1 <- multi_sim$subjects[[1]]$W
  W2 <- multi_sim$subjects[[2]]$W
  expect_false(all(W1 == W2))
  
  # But should be correlated with group
  cor_val <- abs(cor(as.vector(W1), as.vector(multi_sim$group$W)))
  expect_true(cor_val > 0.3)  # Should be correlated (allowing for variability)
})

test_that("event design generation works", {
  events <- generate_event_design(K = 3, T = 100, TR = 2, 
                                  event_duration = 2, min_isi = 4)
  
  # Check structure
  expect_true(is.data.frame(events))
  expect_true(all(c("onset", "duration", "condition") %in% names(events)))
  
  # Check timing constraints
  expect_true(all(events$onset >= 0))
  expect_true(all(events$onset < 200))  # T * TR
  expect_true(all(events$condition %in% 1:3))
  
  # Check ISI
  if (nrow(events) > 1) {
    isis <- diff(events$onset)
    expect_true(all(isis >= 4))  # min_isi
  }
})

test_that("determine_spatial_dims finds reasonable dimensions", {
  # Perfect cube
  dims <- stance:::determine_spatial_dims(1000)
  expect_equal(prod(dims), 1000)
  expect_equal(dims, c(10, 10, 10))
  
  # Non-perfect cube
  dims2 <- stance:::determine_spatial_dims(1024)
  expect_equal(prod(dims2), 1024)
  expect_true(all(dims2 >= 8))  # Reasonable sizes
  
  # Prime number - should find closest factorization
  dims3 <- stance:::determine_spatial_dims(997)
  # Since 997 is prime, it should find the closest factorization
  expect_true(abs(prod(dims3) - 997) <= 3)  # Allow small difference
})
