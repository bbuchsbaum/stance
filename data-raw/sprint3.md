Okay, Sprint 3 will build upon the HRF estimation and Rcpp optimizations from Sprint 2. The main goals for Sprint 3 are:
1.  **Implement Spatial Smoothing for HRFs (GMRF Prior):** This is a key feature for robust HRF estimation.
2.  **Further Rcpp/Performance Optimizations:** Address any remaining R bottlenecks, especially in M-step computations or areas not covered in Sprint 2. Explore some of the advanced speed tips (e.g., FFT for convolution if not fully done, low-rank projection for likelihoods).
3.  **Low-Rank Projection for Likelihoods:** If not already fully implemented as part of S2-T08, make this a primary optimization target.
4.  **Refine User Interface and Diagnostics:** Improve usability and provide more insightful diagnostic tools.
5.  **Prepare for More Complex Simulations/Real Data:** Ensure the model is stable enough for initial tests on more realistic scenarios.

**Sprint 3 Goal:** Implement Gaussian Markov Random Field (GMRF) spatial priors for HRF coefficients to encourage smooth variation across neighboring voxels. Complete the low-rank projection optimization for likelihood calculations. Introduce advanced algorithmic improvements including Stochastic Variational Bayes (SVB) for scalability. Enhance diagnostic visualization leveraging `neuroim2` spatial capabilities. **This sprint builds on Sprint 1a shared infrastructure, Sprint 1 CBD implementation, and Sprint 2 Rcpp optimizations. Maintains focus on CPU-based optimizations per proposal.md performance targets (ROI: 5-15 min, Parcellated: 30-90 min, Full WB: 4-12 hours).** **All code must maintain 0 errors, 0 warnings, ≤1 NOTE on `R CMD check --as-cran` across all platforms.**

**Assumptions for Sprint 3:**
*   **Sprint 1a/1/2 completed:** Shared infrastructure, CBD core, and basic Rcpp optimizations are functional
*   GMRF will use standard 6-connectivity (face-touching neighbors) in 3D voxel space
*   Low-rank projection Y_proj = UᵀY is fully implemented for all likelihood calculations  
*   FFT convolution threshold optimized based on Sprint 2 profiling results
*   SVB is an optional enhancement - core VB must remain stable
*   All spatial operations leverage `neuroim2::NeuroSpace` for coordinate handling
*   **🔹SPEED-TIP #3 from proposal:** Batch voxels for GMRF updates using shared Cholesky factors
*   **CI gate:** GitHub Actions matrix (macOS-latest, ubuntu-latest, windows-latest) must stay green.

---

**Sprint 3 Tickets:**

**Theme 1: Spatially Smoothed HRFs (GMRF Prior)**

1.  **Ticket S3-T01: Implement Spatial Neighborhood Infrastructure** 🟢 **[LEVERAGES NEUROIM2]**
    *   Description: Create spatial connectivity functions leveraging `neuroim2`:
        ```r
        # R/spatial_utils.R
        get_spatial_neighbors <- function(mask, connectivity = 6) {
          # Extract NeuroSpace from mask if NeuroVol
          if (inherits(mask, "NeuroVol")) {
            space <- space(mask)
            coords <- which(mask > 0, arr.ind = TRUE)
          } else {
            # Assume mask is logical array with coordinates
            coords <- which(mask, arr.ind = TRUE)
          }
          
          # Define connectivity offsets
          offsets <- switch(as.character(connectivity),
            "6" = rbind(c(-1,0,0), c(1,0,0), c(0,-1,0), 
                       c(0,1,0), c(0,0,-1), c(0,0,1)),
            "18" = expand.grid(-1:1, -1:1, -1:1)[-14,],  # Remove center
            "26" = expand.grid(-1:1, -1:1, -1:1)[-14,]
          )
          
          # Build adjacency list
          n_voxels <- nrow(coords)
          neighbors <- vector("list", n_voxels)
          
          for (i in 1:n_voxels) {
            # Find valid neighbors
            neighbor_coords <- sweep(offsets, 2, coords[i,], "+")
            # Check bounds and mask
            valid <- apply(neighbor_coords, 1, function(nc) {
              all(nc > 0) && all(nc <= dim(mask)) && 
              mask[nc[1], nc[2], nc[3]] > 0
            })
            
            # Store neighbor indices
            neighbors[[i]] <- which(apply(coords, 1, function(c) {
              any(apply(neighbor_coords[valid,, drop=FALSE], 1, 
                       function(nc) all(nc == c)))
            }))
          }
          
          return(list(neighbors = neighbors, coords = coords))
        }
        
        # Create sparse Laplacian matrix
        create_gmrf_laplacian <- function(neighbors, n_voxels) {
          # Use Matrix package for sparse representation
          i <- j <- x <- integer(0)
          
          for (v in 1:n_voxels) {
            n_neighbors <- length(neighbors[[v]])
            if (n_neighbors > 0) {
              # Diagonal entry
              i <- c(i, v)
              j <- c(j, v) 
              x <- c(x, n_neighbors)
              
              # Off-diagonal entries
              i <- c(i, rep(v, n_neighbors))
              j <- c(j, neighbors[[v]])
              x <- c(x, rep(-1, n_neighbors))
            }
          }
          
          L <- Matrix::sparseMatrix(i = i, j = j, x = x, 
                                   dims = c(n_voxels, n_voxels))
          return(L)
        }
        ```
    *   **🔄 Integrate with neuroim2:** Use `neuroim2::coords()` for voxel-to-world conversions
    *   Add unit tests verifying: Laplacian row sums = 0, symmetry, correct sparsity pattern
    *   Acceptance Criteria: Functions handle both NeuroVol masks and coordinate matrices. Laplacian has correct properties. Performance <1s for 50k voxels. **Full neuroim2 compatibility.**
2.  **Ticket S3-T02: Extend CBD R6 Class for GMRF Support** 🟡 **[EXTENDS SPRINT 1 CBD]**
    *   Description: Update `ContinuousBayesianDecoder` from Sprint 1 to support spatial priors:
        ```r
        # In initialize()
        initialize = function(..., use_gmrf = FALSE, lambda_h = 1.0, 
                            mask = NULL, connectivity = 6) {
          # ... existing Sprint 1 code ...
          
          if (use_gmrf) {
            if (is.null(mask)) {
              stop("Mask required for GMRF prior")
            }
            
            # Build spatial structure using Sprint 3 utilities
            spatial_info <- get_spatial_neighbors(mask, connectivity)
            private$.L_gmrf <- create_gmrf_laplacian(spatial_info$neighbors, 
                                                    length(spatial_info$neighbors))
            private$.lambda_h <- lambda_h
            private$.use_gmrf <- TRUE
            
            # Pre-compute Cholesky factor for efficiency (🔹SPEED-TIP #3)
            # Add small diagonal for numerical stability
            Q <- private$.lambda_h * private$.L_gmrf + 
                 Matrix::Diagonal(nrow(private$.L_gmrf), 1e-6)
            private$.Q_chol <- Matrix::Cholesky(Q)
            
            # Store spatial metadata for visualization
            private$.spatial_info <- spatial_info
          }
        }
        ```
    *   **🔄 Maintain Sprint 1 interface:** Ensure backward compatibility
    *   Store spatial metadata for diagnostic plotting
    *   Acceptance Criteria: CBD initializes with GMRF parameters. Spatial structure correctly extracted from mask. **Extends Sprint 1 implementation cleanly.**
3.  **Ticket S3-T03: Implement GMRF-Regularized HRF Updates**
    *   Description: Modify HRF coefficient updates from Sprint 2 to include spatial smoothing:
        ```r
        private$.update_H_v_gmrf <- function() {
          # Check if using shared utilities from Sprint 1a
          if (!exists("convolve_with_hrf")) {
            source("R/hrf_utils.R")  # Sprint 1a shared utilities
          }
          
          L_basis <- ncol(private$.hrf_basis)
          
          # 🔹SPEED-TIP #3: Batch voxels into blocks with same Laplacian structure
          if (private$.use_batched_gmrf) {
            block_size <- 64  # 4x4x4 blocks
            n_blocks <- ceiling(nrow(private$.H_v) / block_size)
            
            for (block in 1:n_blocks) {
              voxel_idx <- ((block-1)*block_size + 1):min(block*block_size, nrow(private$.H_v))
              
              # Extract block Laplacian (assume same structure within block)
              L_block <- private$.L_gmrf[voxel_idx, voxel_idx]
              Q_block <- Matrix::crossprod(X_block) / private$.sigma2 + 
                        private$.lambda_h * L_block
              
              # Single Cholesky for entire block
              Q_chol_block <- Matrix::Cholesky(Q_block)
              
              # Update all voxels in block
              for (l in 1:L_basis) {
                # ... (rest of update logic)
                private$.H_v[voxel_idx, l] <- Matrix::solve(Q_chol_block, b_block)
              }
            }
          } else {
            # Standard voxel-wise updates (fallback)
            # ... (existing logic)
          }
        }
        ```
    *   **🔄 Use Sprint 1a HRF utilities** for convolution operations
    *   Handle edge cases: singular systems, disconnected voxels
    *   Acceptance Criteria: HRF updates incorporate spatial smoothing. Solver is numerically stable. Smoothness visually apparent. **Batched updates achieve >3x speedup.**
4.  **Ticket S3-T04: Update ELBO for GMRF Prior** 🟡 **[EXTENDS SPRINT 1 ELBO]**
    *   Description: Add spatial prior term to ELBO from Sprint 1:
        ```r
        private$.compute_elbo <- function() {
          # Get base ELBO from Sprint 1 implementation
          elbo_base <- super$.compute_elbo()
          
          if (private$.use_gmrf) {
            # GMRF prior: -0.5 * lambda_h * sum_l (H[:,l]' * L_gmrf * H[:,l])
            gmrf_penalty <- 0
            for (l in 1:ncol(private$.H_v)) {
              h_l <- private$.H_v[,l]
              gmrf_penalty <- gmrf_penalty + 
                as.numeric(t(h_l) %*% private$.L_gmrf %*% h_l)
            }
            
            # Log determinant term (constant if lambda_h fixed)
            # log p(H|lambda) = 0.5 * log|lambda * L| - 0.5 * lambda * H'LH
            # Approximation: log|L| ~ (V-1) * log(lambda) for connected graph
            n_voxels <- nrow(private$.H_v)
            log_det_term <- 0.5 * (n_voxels - 1) * log(private$.lambda_h)
            
            elbo_spatial <- log_det_term - 0.5 * private$.lambda_h * gmrf_penalty
            elbo_base <- elbo_base + elbo_spatial
          }
          
          return(elbo_base)
        }
        ```
    *   Verify ELBO still increases with GMRF
    *   Acceptance Criteria: ELBO correctly includes GMRF prior. Monotonic increase maintained. **Compatible with Sprint 1 convergence checks.**
5.  **Ticket S3-T05: Test GMRF with Spatially Smooth Simulations** 🟢 **[USES SPRINT 1a SIMULATION]**
    *   Description: Extend Sprint 1a simulation framework for spatial validation:
        ```r
        # Extend simulate_fmri_data() from Sprint 1a
        test_gmrf_smoothing <- function() {
          # Generate spatially smooth HRFs
          V <- 1000  # voxels in 10x10x10 cube
          coords <- expand.grid(1:10, 1:10, 1:10)
          
          # Use Sprint 1a simulation with spatial HRFs
          sim_data <- simulate_fmri_data(
            V = V, T = 300, K = 3,
            algorithm = "CBD",  # Use CBD mode from Sprint 1
            spatial_hrf = list(
              correlation_length = 5,  # voxels
              smoothness = 0.8
            )
          )
          
          # Create mask as NeuroVol using neuroim2
          mask_data <- array(TRUE, c(10, 10, 10))
          mask <- neuroim2::NeuroVol(
            data = mask_data,
            space = neuroim2::NeuroSpace(dim = c(10, 10, 10))
          )
          
          # Fit with and without GMRF using Sprint 1 CBD
          cbd_no_gmrf <- ContinuousBayesianDecoder$new(
            Y = sim_data$Y,
            rank = 10,
            n_states = 3,
            hrf_spec = fmrireg::HRF_SPMG1,  # Sprint 1a HRF integration
            use_gmrf = FALSE
          )
          
          cbd_gmrf <- ContinuousBayesianDecoder$new(
            Y = sim_data$Y, 
            rank = 10,
            n_states = 3,
            hrf_spec = fmrireg::HRF_SPMG1,
            use_gmrf = TRUE,
            lambda_h = 10,
            mask = mask
          )
          
          # Fit both models
          cbd_no_gmrf$fit(max_iter = 50)
          cbd_gmrf$fit(max_iter = 50)
          
          # Compare smoothness of estimates
          # GMRF should have lower roughness
          roughness_no_gmrf <- compute_roughness(cbd_no_gmrf$get_hrfs())
          roughness_gmrf <- compute_roughness(cbd_gmrf$get_hrfs())
          
          expect_true(roughness_gmrf < 0.7 * roughness_no_gmrf)
        }
        ```
    *   Visualize HRF maps as NeuroVol objects using `neuroim2::overlay()`
    *   Acceptance Criteria: GMRF reduces roughness by >30%. Visual smoothness apparent. No numerical issues. **Uses Sprint 1a simulation infrastructure.**

**Theme 2: Performance Optimization (Low-Rank Projection & FFT Convolution)**

6.  **Ticket S3-T06: Complete Low-Rank Projection Implementation** 🟡 **[COMPLETES SPRINT 2 OPTIMIZATION]**
    *   Description: Ensure all likelihood calculations use low-rank projection per proposal 🔹SPEED-TIP #2:
        ```cpp
        // In src/likelihood_lowrank.cpp (extends Sprint 2 Rcpp)
        // [[Rcpp::export]]
        arma::mat compute_log_likelihood_lowrank_complete(
            const arma::mat& Y_proj,     // r x T (precomputed U'Y)
            const arma::mat& U,          // V x r (orthonormal)
            const arma::mat& V,          // K x r
            const arma::mat& H_v,        // V x L_basis
            const arma::mat& hrf_basis,  // L_h x L_basis
            const arma::mat& S_gamma,    // K x T (posterior expectations)
            double sigma2) {
          
          int K = V.n_rows;
          int T = Y_proj.n_cols;
          int r = Y_proj.n_rows;
          int V_voxels = U.n_rows;
          
          arma::mat log_lik(K, T);
          double log_const = -0.5 * r * log(2 * M_PI * sigma2);
          
          // Leverage Sprint 1a HRF utilities via Rcpp interface
          // ... (convolution logic)
          
          // 🔹SPEED-TIP #2: Work entirely in r-dimensional space
          #pragma omp parallel for collapse(2)
          for (int t = 0; t < T; t++) {
            for (int k = 0; k < K; k++) {
              // Compute predicted signal in low-rank space
              arma::vec mu_proj_kt = compute_mu_proj(k, t, U, V, H_v, S_conv);
              
              // Squared error in low-rank space
              double sq_error = arma::as_scalar(
                (Y_proj.col(t) - mu_proj_kt).t() * 
                (Y_proj.col(t) - mu_proj_kt)
              );
              
              log_lik(k, t) = log_const - 0.5 * sq_error / sigma2;
            }
          }
          
          return log_lik;
        }
        ```
    *   **Performance target**: 50-100x speedup for V>10k per proposal
    *   Add comprehensive unit tests comparing to full-space calculation
    *   Acceptance Criteria: Low-rank likelihood matches full version within 1e-8. **Achieves target speedup.** Integrates with Sprint 2 Rcpp infrastructure.
7.  **Ticket S3-T07: Optimize FFT Convolution Implementation** 🟡 **[FINALIZES SPRINT 2 FFT WORK]**
    *   Description: Finalize FFT-based convolution strategy per proposal 🔹SPEED-TIP #1:
        ```r
        # R/convolution_utils.R (extends Sprint 1a HRF utilities)
        optimize_convolution_strategy <- function(T, L_h, profile = FALSE) {
          # Empirically determined thresholds based on Sprint 2 profiling
          thresholds <- list(
            fft = 128,      # Use FFT for T > 128 (🔹SPEED-TIP #1)
            batch_fft = 64  # Batch multiple convolutions for T > 64
          )
          
          if (profile) {
            # Run micro-benchmarks to determine optimal threshold
            # Use Sprint 1a convolve_with_hrf() for comparison
            times_direct <- times_fft <- numeric(5)
            
            for (i in 1:5) {
              test_signal <- rnorm(T)
              test_kernel <- fmrireg::HRF_SPMG1(1:min(L_h, 32))
              
              times_direct[i] <- system.time({
                convolve_with_hrf(test_signal, test_kernel, method = "direct")
              })[3]
              
              times_fft[i] <- system.time({
                convolve_with_hrf(test_signal, test_kernel, method = "fft")
              })[3]
            }
            
            # Update threshold based on results
            if (median(times_fft) < median(times_direct)) {
              thresholds$fft <- T
            }
          }
          
          return(thresholds)
        }
        
        # Batch FFT convolution in Rcpp (extends Sprint 2)
        # See src/convolution_fft.cpp
        ```
    *   Cache FFT plans for repeated convolutions
    *   **Performance target**: 10x speedup for T>256 per proposal
    *   Acceptance Criteria: FFT used optimally. **Achieves target speedup.** Correct output for edge cases. Integrates with Sprint 1a HRF framework.
8.  **Ticket S3-T08: Profile and Optimize Remaining Bottlenecks**
    *   Description: Comprehensive profiling of full VB iteration to meet proposal targets:
        ```r
        profile_vb_iteration <- function(cbd_object, target = "parcellated") {
          # Profile against proposal.md performance targets
          targets <- list(
            roi = list(voxels = 1000, time = "5-15 minutes", ram = "<2 GB"),
            parcellated = list(voxels = 5000, time = "30-90 minutes", ram = "<8 GB"),
            full_wb = list(voxels = 50000, time = "4-12 hours", ram = "<16-32GB")
          )
          
          current_target <- targets[[target]]
          
          # Detailed timing of each component
          timings <- profile_cbd_components(cbd_object)
          
          # Check against targets
          total_time <- sum(unlist(timings)) * cbd_object$max_iter
          
          if (total_time > parse_time(current_target$time)) {
            # Identify optimization opportunities
            suggest_optimizations(timings, target)
          }
          
          return(list(
            timings = timings,
            memory = lobstr::mem_used(),
            meets_target = total_time <= parse_time(current_target$time)
          ))
        }
        ```
    *   Create targeted optimization plan for Sprint 4
    *   Acceptance Criteria: Clear performance profile. **Meeting proposal.md targets.** Bottlenecks identified with quantitative metrics.

**Theme 3: Stochastic Variational Bayes (SVB) & Advanced Optimizations**

9.  **Ticket S3-T09: Implement Stochastic Variational Bayes (SVB)** 🟡 **[OPTIONAL - 🔹SPEED-TIP #4]**
    *   Description: Add mini-batch capability for large datasets per proposal:
        ```r
        # Extend Sprint 1 fit() method
        fit = function(max_iter = 100, tol = 1e-4, 
                      batch_size = NULL,  # NULL for full batch
                      learning_rate = function(t) 1/(1 + 0.1*t),
                      ...) {
          
          if (!is.null(batch_size)) {
            # 🔹SPEED-TIP #4: Stochastic VB mode
            T_total <- ncol(private$.Y_data)
            
            # Validate batch size
            if (batch_size < 100 || batch_size > T_total/2) {
              warning("Batch size should be 100-150 for optimal performance")
            }
            
            n_batches <- ceiling(T_total / batch_size)
            
            for (iter in 1:max_iter) {
              # Sample mini-batch (contiguous for temporal coherence)
              batch_start <- sample(1:(T_total - batch_size + 1), 1)
              batch_idx <- batch_start:(batch_start + batch_size - 1)
              
              # E-step on mini-batch using Sprint 1 methods
              Y_batch <- private$.Y_data[, batch_idx]
              Y_proj_batch <- crossprod(private$.U, Y_batch)
              
              log_lik_batch <- private$.compute_log_likelihoods(
                Y_proj_batch, subset = TRUE
              )
              
              # Local forward-backward
              fb_result <- private$.forward_backward(
                log_lik_batch[, batch_idx], 
                private$.Pi, 
                private$.pi0
              )
              
              # Update global parameters with learning rate
              lr <- learning_rate(iter)
              
              # Natural gradient updates
              private$.Pi <- (1 - lr) * private$.Pi + 
                           lr * normalize_rows(fb_result$xi_stats)
              
              # Stochastic updates for U, V using mini-batch statistics
              private$.update_U_V_stochastic(
                Y_batch, fb_result$gamma, lr
              )
              
              # Periodic full evaluation for ELBO
              if (iter %% 10 == 0) {
                elbo <- private$.compute_elbo()
                if (abs(elbo - private$.elbo_prev) < tol) break
              }
            }
          } else {
            # Standard full-batch VB from Sprint 1
            super$fit(max_iter, tol, ...)
          }
        }
        ```
    *   Add variance reduction techniques (momentum, adaptive learning rates)
    *   Acceptance Criteria: SVB converges to similar solution as full VB. **3-5x speedup for T>1000.** Optional feature doesn't break core VB.
10. **Ticket S3-T10: Implement GMRF Batched Solver** 🟡 **[IMPLEMENTS 🔹SPEED-TIP #3]**
    *   Description: Batch voxels for GMRF updates per proposal optimization:
        ```cpp
        // src/gmrf_batch_solver.cpp
        // [[Rcpp::depends(RcppEigen)]]
        #include <RcppEigen.h>
        
        // [[Rcpp::export]]
        Eigen::MatrixXd solve_gmrf_batched(
            const Eigen::SparseMatrix<double>& XtX,
            const Eigen::SparseMatrix<double>& L_gmrf,
            const Eigen::MatrixXd& XtY,  // V x L_basis
            double lambda_h,
            int block_size = 64,  // 🔹SPEED-TIP #3: 8x8x8 blocks
            double tol = 1e-6) {
          
          int V = XtY.rows();
          int L_basis = XtY.cols();
          int n_blocks = (V + block_size - 1) / block_size;
          
          Eigen::MatrixXd H_v(V, L_basis);
          
          // Process voxels in blocks
          #pragma omp parallel for
          for (int b = 0; b < n_blocks; b++) {
            int start = b * block_size;
            int end = std::min(start + block_size, V);
            int block_v = end - start;
            
            // Extract block matrices
            Eigen::SparseMatrix<double> XtX_block = 
              XtX.block(start, start, block_v, block_v);
            Eigen::SparseMatrix<double> L_block = 
              L_gmrf.block(start, start, block_v, block_v);
            
            // Form block precision matrix
            Eigen::SparseMatrix<double> Q_block = 
              XtX_block + lambda_h * L_block;
            
            // Single Cholesky decomposition for entire block
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(Q_block);
            
            // Solve for all basis functions in block
            H_v.block(start, 0, block_v, L_basis) = 
              solver.solve(XtY.block(start, 0, block_v, L_basis));
          }
          
          return H_v;
        }
        ```
    *   Test on problems up to 100k voxels
    *   Acceptance Criteria: **>5x speedup over voxel-wise updates.** Handles ill-conditioned systems gracefully. Memory efficient.

**Theme 4: Diagnostics, Visualization & Integration**

11. **Ticket S3-T11: Enhanced Spatial Visualization** 🟢 **[LEVERAGES NEUROIM2]**
    *   Description: Create comprehensive spatial visualization using `neuroim2`:
        ```r
        # Extend Sprint 1a visualization framework
        plot_spatial_maps <- function(cbd_object, states = NULL, 
                                     slice_coords = NULL,
                                     threshold = 2,
                                     template = NULL) {
          if (is.null(states)) states <- 1:cbd_object$n_states
          
          # Extract spatial maps W_k = U * V[k,]'
          W <- cbd_object$get_spatial_maps()
          
          # Convert to NeuroVol using Sprint 1 metadata
          if (!is.null(cbd_object$private$.neuro_metadata)) {
            space_info <- cbd_object$private$.neuro_metadata$space
            
            plot_list <- list()
            for (k in states) {
              # Create NeuroVol
              vol_k <- neuroim2::NeuroVol(
                data = W[,k],
                space = space_info
              )
              
              # Z-score for visualization
              vol_k_z <- (vol_k - mean(vol_k)) / sd(vol_k)
              
              # Threshold
              vol_k_z[abs(vol_k_z) < threshold] <- 0
              
              # Use neuroim2 overlay functionality
              if (!is.null(template)) {
                p <- neuroim2::overlay(template, vol_k_z,
                                      zlim = c(-4, 4),
                                      alpha = 0.7)
              } else {
                p <- neuroim2::plot(vol_k_z, slices = slice_coords)
              }
              
              plot_list[[k]] <- p
            }
            
            # Combine plots
            do.call(gridExtra::grid.arrange, c(plot_list, ncol = 2))
            
          } else {
            # Fallback to Sprint 1a basic visualization
            plot_spatial_maps_basic(W[,states])
          }
        }
        ```
    *   Support interactive visualization via `neuroim2::view_vol()`
    *   Acceptance Criteria: Publication-quality brain maps. **Full neuroim2 integration.** Works with and without spatial metadata.
12. **Ticket S3-T12: Comprehensive Model Diagnostics** 🟡 **[EXTENDS SPRINT 1]**
    *   Description: Extend Sprint 1 diagnostics with spatial assessments:
        ```r
        # Add to CBD class
        diagnose_model_fit <- function(self, output_dir = NULL) {
          diagnostics <- list(
            # Sprint 1 diagnostics
            convergence = check_convergence_diagnostics(self),
            parameter_recovery = if (exists("true_params")) {
              assess_parameter_recovery(self, true_params)
            },
            
            # Sprint 3 spatial diagnostics
            spatial = if (private$.use_gmrf) list(
              hrf_smoothness = assess_hrf_smoothness(self),
              effective_df = compute_effective_df_gmrf(self),
              spatial_correlation = compute_spatial_autocorrelation(self)
            ),
            
            # Performance diagnostics
            performance = list(
              iteration_time = mean(private$.iteration_times),
              memory_peak = private$.peak_memory,
              meets_targets = check_performance_targets(self)
            )
          )
          
          # Generate HTML report
          if (!is.null(output_dir)) {
            generate_diagnostic_report(diagnostics, output_dir)
          }
          
          return(diagnostics)
        }
        ```
    *   Include residual analysis and goodness-of-fit metrics
    *   Acceptance Criteria: Comprehensive diagnostics available. **Identifies potential issues.** HTML report generation.
13. **Ticket S3-T13: Cross-Validation Framework** 🟡 **[NEW FUNCTIONALITY]**
    *   Description: Implement time-blocked CV for hyperparameter selection:
        ```r
        # R/cross_validation.R
        cbd_cross_validate <- function(Y, mask = NULL,
                                      rank_values = c(5, 10, 15, 20),
                                      lambda_h_values = c(0.1, 1, 10, 100),
                                      n_folds = 5,
                                      metric = "predictive_likelihood",
                                      use_parallel = TRUE) {
          
          # Time-blocked CV (preserve temporal structure)
          T <- ncol(Y)
          fold_size <- floor(T / n_folds)
          fold_starts <- seq(1, T - fold_size + 1, by = fold_size)[1:n_folds]
          
          # Setup parallel backend if requested
          if (use_parallel) {
            future::plan(future::multisession, workers = min(n_folds, 4))
          }
          
          # Grid search using Sprint 1a/1 infrastructure
          results <- expand.grid(
            rank = rank_values,
            lambda_h = lambda_h_values,
            fold = 1:n_folds,
            score = NA
          )
          
          # Parallel CV using future
          cv_results <- future.apply::future_lapply(1:nrow(results), function(i) {
            # ... (CV logic using CBD from Sprint 1)
          })
          
          # Find best parameters
          best_params <- select_best_parameters(cv_results)
          
          return(list(
            results = cv_results,
            best_params = best_params,
            recommendation = recommend_parameters(Y, best_params)
          ))
        }
        ```
    *   Support nested CV for unbiased performance estimation
    *   Acceptance Criteria: CV framework works reliably. Selects reasonable hyperparameters. **Handles edge cases.** Parallelization efficient.

**Theme 5: Integration Testing & Sprint 4 Preparation**

14. **Ticket S3-T14: Comprehensive Integration Testing**
    *   Description: Full system test with all Sprint 1/2/3 features:
        ```r
        test_sprint3_integration <- function() {
          # Use Sprint 1a unified simulation
          set.seed(42)
          sim_data <- simulate_fmri_data(
            V = 5000, T = 600, K = 4,
            algorithm = "CBD",
            spatial_correlation = 5,
            hrf_basis = fmrireg::HRF_BSPLINE,  # Test flexible basis
            snr = 1.5
          )
          
          # Test configurations building on previous sprints
          configs <- list(
            sprint1_baseline = list(
              use_gmrf = FALSE, 
              use_rcpp = FALSE  # Pure R from Sprint 1
            ),
            sprint2_rcpp = list(
              use_gmrf = FALSE, 
              use_rcpp = TRUE   # Sprint 2 optimizations
            ),
            sprint3_gmrf = list(
              use_gmrf = TRUE, 
              use_rcpp = TRUE,
              lambda_h = 10,
              use_batched_gmrf = TRUE  # Sprint 3 batching
            ),
            sprint3_full = list(
              use_gmrf = TRUE, 
              use_rcpp = TRUE,
              lambda_h = 10,
              batch_size = 100,  # SVB mode
              use_lowrank = TRUE,
              use_fft = TRUE
            )
          )
          
          # Verify performance improvements
          results <- run_integration_tests(configs, sim_data)
          
          # Check against proposal targets
          expect_true(results$sprint3_full$time < 
                     parse_minutes("90 minutes"))  # Parcellated target
          
          return(results)
        }
        ```
    *   Verify no memory leaks across platforms
    *   Acceptance Criteria: All configurations work. **Performance gains demonstrated.** No regressions. Meeting proposal targets.
15. **Ticket S3-T15: Documentation Sprint**
    *   Description: Comprehensive documentation update building on Sprint 1a/1:
        - Update package documentation with GMRF features
        - Create vignette: "spatial-priors-with-stance.Rmd"
        - Create vignette: "performance-optimization-guide.Rmd"
        - Document 🔹SPEED-TIPS implementation details
        - Create decision tree for parameter selection
        - Add troubleshooting guide for common issues
    *   **Ensure all examples work with `neuroim2`/`fmrireg` data**
    *   Acceptance Criteria: Documentation complete and clear. **New users can understand and use advanced features.** Examples reproducible.
16. **Ticket S3-T16: Performance Benchmarking Against Proposal Targets**
    *   Description: Verify implementation meets proposal.md performance targets:
        ```r
        benchmark_against_proposal <- function() {
          # Test each scenario from proposal Table
          scenarios <- list(
            roi = list(V = 1000, T = 500, K = 3, r = 10,
                      target = "5-15 minutes", ram = "2GB"),
            parcellated = list(V = 5000, T = 500, K = 4, r = 15,
                             target = "30-90 minutes", ram = "8GB"),
            full_wb = list(V = 50000, T = 1000, K = 5, r = 20,
                          target = "4-12 hours", ram = "32GB")
          )
          
          results <- list()
          for (name in names(scenarios)) {
            sc <- scenarios[[name]]
            
            # Generate appropriate test data
            Y <- simulate_fmri_data(V = sc$V, T = sc$T, K = sc$K,
                                   algorithm = "CBD")$Y
            
            # Create optimized CBD instance
            cbd <- ContinuousBayesianDecoder$new(
              Y = Y,
              rank = sc$r,
              n_states = sc$K,
              use_gmrf = (name != "roi"),
              use_rcpp = TRUE,
              optimization_level = "aggressive"
            )
            
            # Time the fit
            time_taken <- system.time({
              cbd$fit(max_iter = 50)
            })["elapsed"]
            
            # Check memory
            mem_used <- lobstr::mem_used()
            
            results[[name]] <- list(
              time_minutes = time_taken / 60,
              memory_gb = mem_used / 1e9,
              meets_target = time_taken < parse_time_minutes(sc$target),
              efficiency = sc$V * sc$T / time_taken  # voxel-TRs per second
            )
          }
          
          return(results)
        }
        ```
    *   Document any deviations from targets with explanations
    *   Acceptance Criteria: **Clear performance measurements against proposal targets.** Explanations for any gaps.
17. **Ticket S3-T17: Sprint 4 Planning & Risk Assessment**
    *   Description: Plan for real data application in Sprint 4:
        - Review all Sprint 1/2/3 components for production readiness
        - Identify remaining gaps for real fMRI analysis
        - Create risk mitigation plan for Sprint 4
        - Draft real-data preprocessing pipeline
        - Plan collaboration with neuroimaging labs for validation
    *   **Consider integration with fMRIPrep outputs via `neuroim2`**
    *   Acceptance Criteria: Sprint 3 code reviewed. **Clear Sprint 4 roadmap.** Risk mitigation strategies defined.

---

**Cross-Cutting Risks & Mitigations:**

| Risk | Impact | Mitigation |
|------|--------|------------|
| GMRF solver instability with disconnected regions | High | Add connected component detection, solve separately |
| Memory explosion with 50k+ voxel GMRF | High | Implement out-of-core solving, use approximate methods |
| SVB convergence issues | Medium | Make optional, provide convergence diagnostics |
| Platform-specific sparse solver failures | Medium | Provide multiple solver backends, test extensively |
| Performance regression with new features | High | Continuous benchmarking, feature flags for optimizations |
| Batch size selection for SVB | Low | Auto-tuning based on T and available memory |

---

**Exit Criteria for Sprint 3:**
1.  **Functionality**: All spatial prior and optimization features working
    - GMRF implementation reduces HRF roughness by >30%
    - Low-rank likelihood computation matches full version (diff <1e-8)
    - FFT convolution provides 10x speedup for T>256
    - SVB converges to similar solution as full VB (optional)
2.  **Performance**: Meeting computational targets from proposal.md
    - ROI (1k voxels): 5-15 minutes ✓
    - Parcellated (5k voxels): 30-90 minutes ✓
    - Full WB (50k voxels): 4-12 hours (overnight) ✓
    - Memory usage within specified bounds
3.  **Quality**: Maintaining software standards
    - R CMD check: 0 errors, 0 warnings, ≤1 NOTE
    - Test coverage remains >85%
    - All CI builds green including sparse solver tests
    - No performance regressions from Sprint 2
4.  **Integration**: Building on previous sprints
    - Seamless integration with Sprint 1a shared infrastructure
    - Sprint 1 CBD implementation properly extended
    - Sprint 2 optimizations preserved and enhanced
    - Full `neuroim2`/`fmrireg` ecosystem compatibility
5.  **Documentation**: Ready for advanced users and Sprint 4
    - Spatial prior vignette with examples
    - Performance optimization guide with benchmarks
    - Clear guidance on parameter selection
    - Sprint 4 plan addressing real data challenges

**Sprint 3 Summary:** This sprint delivers crucial spatial regularization and performance optimizations that make the CBD implementation practical for real neuroimaging applications. By implementing GMRF priors, completing low-rank optimizations, and adding optional SVB, we achieve the computational targets outlined in the proposal while maintaining the scientific rigor and software quality established in previous sprints. The implementation remains CPU-focused and builds systematically on the foundation from Sprints 1a, 1, and 2.