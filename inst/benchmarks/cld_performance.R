# Performance Benchmark Suite for CLD (S1a-T16)

library(stance)
library(microbenchmark)
library(ggplot2)

# Function to run benchmarks
run_cld_benchmarks <- function() {
  cat("CLD Performance Benchmark Suite\n")
  cat("===============================\n\n")
  
  # Test configurations
  configs <- list(
    small = list(V = 1000, T = 200, K = 3, name = "Small (ROI)"),
    medium = list(V = 5000, T = 500, K = 5, name = "Medium (Parcellated)"),
    large = list(V = 20000, T = 300, K = 4, name = "Large (Whole Brain)")
  )
  
  results <- list()
  
  for (config_name in names(configs)) {
    config <- configs[[config_name]]
    cat(sprintf("\nBenchmarking %s: V=%d, T=%d, K=%d\n", 
                config$name, config$V, config$T, config$K))
    cat(rep("-", 50), "\n", sep = "")
    
    # Generate test data
    set.seed(123)
    sim <- simulate_fmri_data(
      V = config$V, 
      T = config$T, 
      K = config$K,
      algorithm = "CLD",
      verbose = FALSE
    )
    
    # Benchmark key operations
    
    # 1. Data validation
    time_validate <- microbenchmark(
      validate = validate_fmri_input(sim$Y, verbose = FALSE),
      times = 10
    )
    
    # 2. HRF setup and convolution
    hrf <- setup_hrf_kernel("spmg1")
    time_convolution <- microbenchmark(
      convolution = convolve_with_hrf(sim$S, hrf),
      times = 10
    )
    
    # 3. CLD initialization (includes W learning)
    time_init <- microbenchmark(
      initialization = ContinuousLinearDecoder$new(
        Y = sim$Y,
        S_design = sim$S,
        verbose = FALSE
      ),
      times = 5
    )
    
    # 4. Memory usage
    cld <- ContinuousLinearDecoder$new(
      Y = sim$Y,
      S_design = sim$S,
      verbose = FALSE
    )
    
    if (requireNamespace("pryr", quietly = TRUE)) {
      mem_usage <- pryr::object_size(cld)
    } else {
      mem_usage <- object.size(cld)
    }
    
    # Store results
    results[[config_name]] <- list(
      config = config,
      validate = summary(time_validate),
      convolution = summary(time_convolution),
      initialization = summary(time_init),
      memory = mem_usage
    )
    
    # Print results
    cat("\nOperation timings (median):\n")
    cat(sprintf("  Data validation: %.3f ms\n", 
                median(time_validate$time) / 1e6))
    cat(sprintf("  HRF convolution: %.3f ms\n", 
                median(time_convolution$time) / 1e6))
    cat(sprintf("  CLD initialization: %.3f s\n", 
                median(time_init$time) / 1e9))
    cat(sprintf("  Memory usage: %.1f MB\n", 
                as.numeric(mem_usage) / 1024^2))
  }
  
  # Test FISTA performance (if compiled)
  if (exists("fista_tv_rcpp")) {
    cat("\n\nFISTA Optimization Benchmarks\n")
    cat(rep("=", 30), "\n\n", sep = "")
    
    # Small problem for FISTA testing
    sim_fista <- simulate_fmri_data(V = 1000, T = 100, K = 3, verbose = FALSE)
    cld_fista <- ContinuousLinearDecoder$new(
      Y = sim_fista$Y,
      S_design = sim_fista$S,
      verbose = FALSE
    )
    
    time_fit <- microbenchmark(
      fista_fit = cld_fista$fit(max_iter = 20, verbose = FALSE),
      times = 3
    )
    
    cat(sprintf("FISTA fit (20 iterations): %.3f s\n", 
                median(time_fit$time) / 1e9))
  }
  
  # Compare matrix vs neuroim2 input
  cat("\n\nData Structure Comparison\n")
  cat(rep("=", 25), "\n\n", sep = "")
  
  # Standard matrix
  time_matrix <- microbenchmark(
    matrix_input = {
      cld <- ContinuousLinearDecoder$new(
        Y = sim$Y,
        S_design = sim$S,
        verbose = FALSE
      )
    },
    times = 5
  )
  
  # Create mock neuroim2 structure
  Y_neuroim <- structure(
    sim$Y,
    class = c("NeuroVec", "array"),
    space = list(dim = c(10, 10, 10, config$T))
  )
  
  # Time with neuroim2 input (would work with real neuroim2)
  # time_neuroim <- microbenchmark(
  #   neuroim_input = {
  #     cld <- ContinuousLinearDecoder$new(
  #       Y = Y_neuroim,
  #       S_design = sim$S,
  #       verbose = FALSE
  #     )
  #   },
  #   times = 5
  # )
  
  cat(sprintf("Matrix input: %.3f s\n", median(time_matrix$time) / 1e9))
  # cat(sprintf("NeuroVec input: %.3f s\n", median(time_neuroim$time) / 1e9))
  
  # Create performance plots
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Scaling plot
    scaling_data <- data.frame(
      V = sapply(results, function(x) x$config$V),
      T = sapply(results, function(x) x$config$T),
      K = sapply(results, function(x) x$config$K),
      init_time = sapply(results, function(x) median(x$initialization$time)) / 1e9,
      memory_mb = sapply(results, function(x) as.numeric(x$memory)) / 1024^2
    )
    
    p1 <- ggplot(scaling_data, aes(x = V, y = init_time)) +
      geom_point(size = 3) +
      geom_line() +
      scale_x_log10() +
      labs(x = "Number of Voxels", 
           y = "Initialization Time (s)",
           title = "CLD Initialization Time Scaling") +
      theme_minimal()
    
    p2 <- ggplot(scaling_data, aes(x = V, y = memory_mb)) +
      geom_point(size = 3) +
      geom_line() +
      scale_x_log10() +
      labs(x = "Number of Voxels", 
           y = "Memory Usage (MB)",
           title = "CLD Memory Usage Scaling") +
      theme_minimal()
    
    # Save plots
    ggsave("cld_time_scaling.png", p1, width = 6, height = 4)
    ggsave("cld_memory_scaling.png", p2, width = 6, height = 4)
    
    cat("\nPlots saved: cld_time_scaling.png, cld_memory_scaling.png\n")
  }
  
  return(results)
}

# Function to profile specific operations
profile_cld_operations <- function(V = 5000, T = 200, K = 4) {
  cat("\nDetailed Operation Profiling\n")
  cat("============================\n\n")
  
  # Generate data
  sim <- simulate_fmri_data(V = V, T = T, K = K, verbose = FALSE)
  hrf <- setup_hrf_kernel("spmg1")
  
  # Profile individual steps
  
  # 1. Convolution methods
  cat("Convolution method comparison:\n")
  
  # Direct convolution (small T)
  X_small <- sim$S[, 1:50]
  time_direct <- microbenchmark(
    direct = convolve_with_hrf(X_small, hrf, method = "direct"),
    times = 20
  )
  
  # FFT convolution (large T)
  time_fft <- microbenchmark(
    fft = convolve_with_hrf(sim$S, hrf, method = "fft"),
    times = 20
  )
  
  cat(sprintf("  Direct (T=50): %.3f ms\n", median(time_direct$time) / 1e6))
  cat(sprintf("  FFT (T=%d): %.3f ms\n", T, median(time_fft$time) / 1e6))
  
  # 2. Block processing in GLM
  if (exists("speedglm::speedlm.fit")) {
    cat("\nGLM fitting methods:\n")
    
    # Prepare convolved design
    X_conv <- t(convolve_with_hrf(sim$S, hrf))
    
    # Standard lm.fit
    time_lm <- microbenchmark(
      lm.fit = {
        for (v in 1:min(100, V)) {
          lm.fit(x = X_conv, y = sim$Y[v, ])
        }
      },
      times = 5
    )
    
    # speedglm
    time_speedglm <- microbenchmark(
      speedglm = {
        for (v in 1:min(100, V)) {
          speedglm::speedlm.fit(X = X_conv, y = sim$Y[v, ], intercept = FALSE)
        }
      },
      times = 5
    )
    
    cat(sprintf("  lm.fit (100 voxels): %.3f ms\n", 
                median(time_lm$time) / 1e6))
    cat(sprintf("  speedglm (100 voxels): %.3f ms\n", 
                median(time_speedglm$time) / 1e6))
  }
  
  # 3. SVD performance
  cat("\nSVD performance:\n")
  B <- matrix(rnorm(V * K), V, K)
  
  time_svd_full <- microbenchmark(
    full_svd = svd(B),
    times = 5
  )
  
  time_svd_truncated <- microbenchmark(
    truncated_svd = svd(B, nu = min(20, K), nv = min(20, K)),
    times = 5
  )
  
  cat(sprintf("  Full SVD: %.3f ms\n", median(time_svd_full$time) / 1e6))
  cat(sprintf("  Truncated SVD (r=20): %.3f ms\n", 
              median(time_svd_truncated$time) / 1e6))
}

# Function to test parallelization benefits
test_parallel_scaling <- function() {
  if (!requireNamespace("future", quietly = TRUE)) {
    cat("Package 'future' not available for parallel testing\n")
    return(NULL)
  }
  
  cat("\nParallel Scaling Test\n")
  cat("====================\n\n")
  
  library(future)
  library(future.apply)
  
  # Test data
  V <- 10000
  n_blocks <- 10
  block_size <- V / n_blocks
  
  Y <- matrix(rnorm(V * 100), V, 100)
  X <- matrix(rnorm(3 * 100), 3, 100)
  
  # Sequential processing
  plan(sequential)
  time_seq <- system.time({
    results_seq <- future_lapply(1:n_blocks, function(b) {
      start_idx <- (b - 1) * block_size + 1
      end_idx <- b * block_size
      
      # Simulate GLM fitting
      block_results <- matrix(0, block_size, 3)
      for (v in 1:block_size) {
        block_results[v, ] <- rnorm(3)  # Simulate coefficients
      }
      block_results
    })
  })
  
  # Parallel processing
  plan(multisession, workers = 4)
  time_par <- system.time({
    results_par <- future_lapply(1:n_blocks, function(b) {
      start_idx <- (b - 1) * block_size + 1
      end_idx <- b * block_size
      
      # Simulate GLM fitting
      block_results <- matrix(0, block_size, 3)
      for (v in 1:block_size) {
        block_results[v, ] <- rnorm(3)  # Simulate coefficients
      }
      block_results
    })
  })
  
  plan(sequential)  # Reset
  
  cat(sprintf("Sequential: %.3f s\n", time_seq["elapsed"]))
  cat(sprintf("Parallel (4 cores): %.3f s\n", time_par["elapsed"]))
  cat(sprintf("Speedup: %.2fx\n", time_seq["elapsed"] / time_par["elapsed"]))
}

# Run all benchmarks
if (interactive()) {
  results <- run_cld_benchmarks()
  profile_cld_operations()
  test_parallel_scaling()
}