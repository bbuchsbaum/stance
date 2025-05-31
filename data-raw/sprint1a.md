**Sprint 1a (Revised): Continuous Linear Decoder (CLD) Prototype**

**Sprint 1a Goal:** Implement a functional R6 class for the Continuous Linear Decoder (CLD) capable of learning spatial maps (`W`) via GLM+SVD (using user-provided `S_design`) and estimating "soft" state activations (`X_hat`) via FISTA with a robust Total Variation (TV) penalty (Condat's algorithm). The implementation will use R for high-level logic and RcppArmadillo for the FISTA loop and gradient computation. Build on established R neuroimaging infrastructure (`neuroim2`, `fmrireg`) for data handling and HRF components. **This sprint serves as the foundation for Sprint 1 (CBD), establishing ~60-70% of the shared infrastructure including project setup, neuroimaging data integration, HRF utilities, simulation framework, and testing infrastructure.**

**Assumptions for Sprint 1a:**
*   Input `Y` is a V x T matrix (voxels x time) or `neuroim2::NeuroVec` object.
*   A single, shared canonical HRF is used, leveraging `fmrireg` HRF basis functions.
*   User provides a `K x T` design matrix `S_design` for learning `W`.
*   Integration with `neuroim2` data structures for neuroimaging data compatibility.
*   **Foundation for Sprint 1:** All shared infrastructure components will be designed for reuse in the Continuous Bayesian Decoder implementation.

---

**Sprint 1a (Revised): Continuous Linear Decoder (CLD) - Granular Tickets**

**Theme 1: Shared Infrastructure Setup (🟢 REUSABLE FOR SPRINT 1)**

1.  **Ticket S1a-T01 (Revised): Initialize Project Repository & Shared Environment** 🟢 **[SHARED FOUNDATION]**
    *   Description: Set up Git repository. Define R project structure (R package layout). Establish `renv`. **Add `usethis::use_github_action("check-standard")` for automated R CMD checks.**
        *   Include `neuroim2` and `fmrireg` as key dependencies in the project setup.
        *   Configure package structure to support neuroimaging data workflows.
        *   **🔄 Create shared directories:** `R/`, `src/`, `tests/testthat/`, `inst/`
        *   **🔄 Establish shared files:** `R/stance-package.R`, `R/data_structures.R`, `R/hrf_utils.R`, `R/simulate.R`, `R/visualization.R`
    *   Acceptance Criteria: Project initialized, `renv.lock` created, GitHub Action for checks configured. Neuroimaging dependencies properly declared. **Shared infrastructure files created.**
2.  **Ticket S1a-T02 (Revised): Implement Shared Data Structure Utilities** 🟢 **[SHARED FOUNDATION]**
    *   Description: Create `R/data_structures.R` with functions for seamless matrix ↔ `neuroim2::NeuroVec`/`neuroim2::NeuroVol` conversion:
        *   `validate_fmri_input()`: Check dimensions, detect data type, validate spatial metadata
        *   `extract_data_matrix()`: Convert neuroimaging objects to matrices for computation
        *   `restore_spatial_structure()`: Convert results back to neuroimaging objects with preserved metadata
        *   `check_temporal_alignment()`: Validate TR alignment between data and design
        *   **🔄 Design for both CLD and CBD:** Ensure functions work for both deterministic and probabilistic outputs
    *   Acceptance Criteria: Robust data structure handling for both algorithms. Compatible with `neuroim2` ecosystem. **Shared utilities ready for Sprint 1 reuse.**
3.  **Ticket S1a-T03 (Revised): Implement Shared HRF Utilities** 🟢 **[SHARED FOUNDATION]**
    *   Description: Create `R/hrf_utils.R` with HRF processing functions leveraging `fmrireg`:
        *   `setup_hrf_kernel()`: Accept `fmrireg` HRF specs or numeric vectors, normalize to unit energy
        *   `convolve_with_hrf()`: Efficient convolution with FFT fallback thresholds
        *   `hrf_basis_matrix()`: Generate basis matrices for flexible HRF modeling
        *   `validate_hrf_spec()`: Check HRF specifications and provide informative errors
        *   **🔄 Design for CBD HRF modeling:** Support both fixed HRFs (CLD) and voxel-wise HRFs (CBD)
    *   Acceptance Criteria: HRF utilities compatible with both algorithms. Integration with `fmrireg` ecosystem. **Ready for Sprint 1 HMM convolution operations.**
4.  **Ticket S1a-T04 (Revised): Implement Unified Simulation Framework** 🟢 **[SHARED FOUNDATION]**
    *   Description: Create `R/simulate.R` with comprehensive simulation supporting both algorithms:
        *   `simulate_fmri_data()`: Unified function with algorithm-specific options
        *   Support CLD mode: GLM-based spatial maps, continuous activations with TV structure
        *   Support CBD mode: Low-rank spatial maps, Markov chain state sequences
        *   **🔄 Shared parameters:** `V`, `T`, `K`, `rank`, `hrf_spec`, `snr`, neuroimaging output options
        *   Return format compatible with both CLD and CBD workflows
        *   Include realistic neuroimaging characteristics (spatial coherence, temporal structure)
    *   Acceptance Criteria: Single simulation framework supporting both algorithms. **Enables direct CLD vs CBD comparison in Sprint 1.**

**Theme 2: CLD-Specific R6 Class Structure**

5.  **Ticket S1a-T05 (Revised): Define `ContinuousLinearDecoder` R6 Class Skeleton**
    *   Description: Create `ContinuousLinearDecoder.R`. Define R6 class with private fields (e.g., `$.Y`, `$.hrf`, `$.rank`, `$.lambda_tv`, `$.W`, `$.WtY`, `$.X_hat`, `$.config`, `$.L_fista`) and public methods. **Create `R/cld-math.R` for small math helper functions.**
        *   Design class to handle both standard matrices and `neuroim2::NeuroVec`/`neuroim2::NeuroVol` objects seamlessly using shared utilities.
        *   Include fields for `neuroim2` spatial metadata when working with neuroimaging data.
        *   **🔄 Design patterns for CBD reuse:** Use similar R6 structure and method signatures that CBD can follow
    *   Acceptance Criteria: R6 class definition exists. Helper file `cld-math.R` created. Compatible with neuroimaging data structures. **Provides R6 template for Sprint 1 CBD class.**
6.  **Ticket S1a-T06 (Revised): Implement R6 `initialize` Method (Core Parameters & `W` Learning)**
    *   Description: Flesh out the `initialize` method using shared utilities:
        *   Inputs: `Y` (matrix or `neuroim2::NeuroVec`), `S_design` (`K x T` matrix), `hrf_kernel` (or `fmrireg` HRF specification), `rank`, `lambda_tv`.
        *   Use `validate_fmri_input()` and `setup_hrf_kernel()` from shared utilities
        *   Validate `S_design` (dims, values; warn if `colMeans(S_design)` suggest missing convolution by user for event designs).
        *   **Auto-detect and handle `neuroim2` data structures** using shared conversion functions
        *   Store `Y`, `S_design`, normalized `hrf_kernel`, `rank`, `lambda_tv`.
        *   Call `private$.learn_W_via_glm_svd(self$.Y, self$.S_design, self$.hrf, self$.rank)` to compute and store `self$.W`.
        *   **Pre-compute and cache `self$.WtY = crossprod(self$.W, self$.Y)`.**
        *   Initialize `self$.X_hat` (e.g., to zeros `K x T`).
        *   (Placeholder) Estimate Lipschitz constant `L` for FISTA and store in `private$.L_fista`.
    *   Acceptance Criteria: CLD object instantiated using shared infrastructure. `W` learned, `WtY` cached. `S_design` validated. `X_hat` initialized. **Demonstrates initialization patterns for Sprint 1 CBD.**
7.  **Ticket S1a-T07 (Revised): Implement R6 Getter Methods & Shared S3 Interface** 🟢 **[PARTIALLY SHARED]**
    *   Description: Implement getter methods and start shared visualization framework:
        *   `get_spatial_maps()`: Returns `self$.W` with optional `neuroim2::NeuroVol` conversion using shared utilities
        *   `get_activations()`: Returns `self$.X_hat`
        *   **🔄 Start `R/visualization.R`** with shared plotting functions for both algorithms:
            - `plot_spatial_maps()`: Brain map visualization using `neuroim2` integration
            - `plot_convergence()`: Generic convergence plotting (FISTA objective / ELBO)
            - `plot_state_timecourse()`: State activation over time
        *   **Implement `as.list.ContinuousLinearDecoder` S3 method** for inspection
    *   Acceptance Criteria: CLD getters functional. **Shared visualization framework started for Sprint 1 reuse.**

**Theme 3: CLD-Specific Algorithms (FISTA + TV)**

8.  **Ticket S1a-T08 (Revised): Implement `private$.learn_W_via_glm_svd` (Optimized R)**
    *   Description: This private method takes `Y_input`, `S_design_input`, `hrf_kernel_input`, and `rank_input`.
        1.  **Create Convolved Design Matrix:** Convolve each row of `S_design_input` (K regressors) with `hrf_kernel_input` (using `fmrireg` convolution functions where applicable) to get `X_conv` (`K x T`).
        2.  **Run Voxel-wise GLMs:** For each voxel, run a GLM: `Y_v,: ~ t(X_conv)[,1] + ... + t(X_conv)[,K]`. Collect `K` beta coefficients per voxel into `B` (`V x K`). **Use `speedglm::speedglm.wfit` for GLMs. If `V > 5000`, process voxels in chunks (e.g., using `data.table` or `lapply` over blocks).**
        3.  **Truncated SVD:** Perform SVD: `B_svd <- svd(B, nu = rank_input, nv = rank_input)` for efficiency.
        4.  **Construct rank-r `W`:** 
            *   Extract: `U_r = B_svd$u[, 1:rank_input]` (V x r), `S_r = B_svd$d[1:rank_input]` (r values), `V_r = B_svd$v[, 1:rank_input]` (K x r).
            *   **Construct rank-r approximation: `W_learned = U_r %*% diag(S_r) %*% t(V_r)`** which gives a V x K matrix of rank r.
            *   Note: This is the low-rank approximation of B, which captures the main spatial patterns.
        *   Handle both standard matrix and `neuroim2::NeuroVec` inputs efficiently, maintaining spatial information when available.
    *   Acceptance Criteria: `self$.W` (a `V x K` matrix of effective rank `r`) is computed and stored. Uses `speedglm` and handles large `V`. Works with neuroimaging data structures and `fmrireg` HRF components.
    *   *Note:* For unsupervised starts if `S_design` is not provided, default to a KxT identity pulse train (K=rank) as `S_design` and document this.
9.  **Ticket S1a-T09 (Revised): Implement FFT-based Convolution Operator (R, with `stats::convolve` fallback)**
    *   Description: Create `convolve_rows_with_hrf(matrix_X, hrf_kernel, use_fft_threshold = 256)`. **If `ncol(matrix_X) > use_fft_threshold`, use `stats::mvfft` for efficient FFT-based convolution. Otherwise, use `apply(matrix_X, 1, function(row) stats::convolve(row, hrf_kernel, conj = FALSE, type = "filter")[1:length(row)])` for direct convolution.** Ensure output dimensions are correct (`K x T`).
        *   Consider integration with `fmrireg` convolution functions for consistency with neuroimaging conventions.
        *   Ensure compatibility when `hrf_kernel` is specified via `fmrireg` HRF basis functions.
    *   Acceptance Criteria: Function performs row-wise convolution, choosing FFT or direct method based on `T`. Compatible with `fmrireg` HRF specifications.
10. **Ticket S1a-T10: Implement Gradient for FISTA (RcppArmadillo)**
    *   Description: RcppArmadillo function `compute_gradient_fista_rcpp(Y_or_WtY_arg, W_arg, H_star_X_arg, hrf_kernel_arg, precomputed_WtY = FALSE)`.
        *   If `precomputed_WtY` is TRUE, `Y_or_WtY_arg` is `WtY` (K x T).
        *   `H_star_X_arg` is `K x T` matrix `convolve_rows_with_hrf(X_current, hrf_kernel)`.
        *   If `precomputed_WtY` is FALSE:
            - `Residual = Y_arg - W_arg * H_star_X_arg` (V x T)
            - `Grad_term = crossprod(W_arg, Residual)` (K x T)
        *   If `precomputed_WtY` is TRUE:
            - `WtW_H_star_X = crossprod(W_arg) * H_star_X_arg` (K x T)  
            - `Grad_term = Y_or_WtY_arg - WtW_H_star_X` (K x T)
        *   `Grad_L2 = -convolve_rows_with_hrf_transposed(Grad_term, hrf_kernel_arg)` (where `_transposed` uses time-reversed HRF).
        *   **Optimize with in-place Armadillo operations. Pre-compute and cache `crossprod(W_arg)` if used repeatedly.**
        *   Ensure efficient handling when input data originates from `neuroim2` structures (converted to matrices for Rcpp processing).
    *   Acceptance Criteria: Rcpp function computes gradient efficiently, leveraging `WtY` cache. Compatible with neuroimaging data workflows.
11. **Ticket S1a-T11 (Revised): Implement Proximal Operator for TV Penalty (Condat's Algorithm in Rcpp)**
    *   Description: Write an Rcpp function `prox_tv_condat_rcpp(X_tilde_row, lambda_tv_step)` that implements Condat's direct algorithm for the proximal operator of the Total Variation penalty. **This computes `argmin_z {0.5 ||z - X_tilde_row||^2 + lambda_tv_step * TV(z)}` where TV(z) = sum(|z[i] - z[i-1]|).** Source a permissively licensed C/C++ implementation of Condat's algorithm and integrate into Rcpp. This operates row-wise on `X_tilde`.
    *   Acceptance Criteria: Rcpp function applies robust 1D TV proximal step for each row of `X_tilde`. Unit test with ramp+spike input. Performance optimized for neuroimaging time series lengths.
12. **Ticket S1a-T12 (Revised): Implement `private$.fista_tv` Loop (RcppArmadillo)**
    *   Description: Implement FISTA in RcppArmadillo.
        *   **Lipschitz constant estimation:** 
            - Option 1 (simple): `L = ||W||_F^2` where ||·||_F is Frobenius norm, assuming HRF is normalized.
            - Option 2 (tighter): Use power iteration to estimate the spectral norm of the linear operator A where A(X) = H^T * (W^T * W * (H * X)).
            - Store chosen estimate in `private$.L_fista`.
        *   Use `self$.WtY` to speed up gradient calculation.
        *   Call `prox_tv_condat_rcpp`.
        *   Implement standard FISTA momentum updates: `t_{k+1} = (1 + sqrt(1 + 4*t_k^2))/2`, `Z_{k+1} = X_{k+1} + ((t_k - 1)/t_{k+1}) * (X_{k+1} - X_k)`.
        *   Include convergence check based on relative change in objective or X.
        *   Ensure efficient memory usage when working with neuroimaging-sized datasets.
    *   Acceptance Criteria: FISTA loop in Rcpp uses Condat's TV prox and precomputed `L`, `WtY`. Optimized for neuroimaging data scales.
13. **Ticket S1a-T13 (Revised): Implement R6 `fit` Method (Calls FISTA)** (Enhanced for neuroimaging compatibility)
    *   Description: Implement the public `fit` method that calls the private FISTA solver.
        *   Ensure method works seamlessly regardless of whether input was standard matrix or `neuroim2` data structure.
        *   Provide progress reporting suitable for neuroimaging analysis workflows (e.g., every 10 iterations).
        *   Return self invisibly for method chaining: `invisible(self)`.
    *   Acceptance Criteria: Public interface for FISTA fitting with neuroimaging workflow compatibility.

**Theme 4: Shared Testing & Validation Infrastructure** 🟢 **[SHARED FOUNDATION]**

14. **Ticket S1a-T14 (Revised): Comprehensive Shared Test Suite** 🟢 **[SHARED FOUNDATION]**
    *   Description: Create test infrastructure for both algorithms:
        *   `tests/testthat/test-infrastructure.R`: Test shared utilities (data structures, HRF, simulation)
        *   `tests/testthat/test-cld.R`: CLD-specific algorithm tests
        *   **🔄 Design test fixtures** for Sprint 1 reuse: small datasets, known solutions, edge cases
        *   **🔄 Create test utilities** for parameter recovery metrics, spatial correlation, temporal accuracy
        *   Test neuroimaging data structure preservation throughout pipelines
        *   Include performance benchmarks for both matrix and `neuroim2` workflows
    *   Acceptance Criteria: Robust shared testing framework. **Ready for Sprint 1 CBD test integration.**
15. **Ticket S1a-T15 (Revised): Documentation & Package Polish** 🟢 **[SHARED FOUNDATION]**
    *   Description: Complete documentation emphasizing shared infrastructure:
        *   Roxygen2 documentation for all shared utilities and CLD-specific functions
        *   **🔄 Create shared vignette framework:** `vignettes/getting-started.Rmd` structure for both algorithms
        *   Document neuroimaging integration patterns for Sprint 1 reuse
        *   **🔄 Package structure:** Clean NAMESPACE, proper imports, shared dependencies
        *   **🔄 Add infrastructure for Sprint 1:** Document how CBD will extend CLD foundation
    *   Acceptance Criteria: Well-documented shared infrastructure. **Clear integration points for Sprint 1 CBD development.**

**Theme 5: Low-Effort Adds & Risk Mitigation**

16. **Ticket S1a-T16: Add Performance Benchmark Suite**
    *   Description: Create benchmark script `inst/benchmarks/cld_performance.R`:
        *   Use `microbenchmark` to time key operations (W learning, FISTA iterations).
        *   Test across different problem sizes: (V,T,K) = {(1000,200,3), (5000,500,5), (20000,300,4)}.
        *   Compare performance with standard matrices vs `neuroim2` objects.
        *   Profile memory usage with `profmem::profmem()`.
        *   Generate performance report with timings and scaling behavior.
        *   Benchmark performance with both standard matrices and `neuroim2` data structures.
        *   Include timing comparisons for different data input types and sizes typical of neuroimaging.
    *   Acceptance Criteria: Basic benchmark script available for both data structure types and neuroimaging-relevant scales.
17. **Ticket S1a-T17: Add Visualization Methods**
    *   Description: Implement visualization methods:
        *   `plot.ContinuousLinearDecoder()`: Creates multi-panel plot showing predicted states over time.
        *   `plot_spatial_maps()`: Visualizes learned W as brain maps when neuroimaging data is used.
        *   `plot_convergence()`: Shows FISTA objective function decrease over iterations.
        *   Integration with `neuroim2::overlay()` for spatial map visualization.
        *   When spatial information is available, provide basic spatial map visualization using `neuroim2` capabilities.
        *   Include HRF visualization using `fmrireg` plotting conventions where applicable.
    *   Acceptance Criteria: Simple plot method available with neuroimaging-appropriate visualizations.
18. **Ticket S1a-T18: Robustness Testing and Edge Cases**
    *   Description: Create comprehensive test suite for edge cases:
        *   **Condat TV prox tests:** 
            - Constant signal (should be unchanged)
            - Linear ramp (should be unchanged for small lambda)
            - Step function (should be smoothed)
            - Impulse noise on smooth signal (should be removed)
        *   **Numerical stability tests:**
            - Very high/low SNR data
            - Poorly conditioned W (high condition number)
            - Long time series (T > 10000)
        *   **Data structure tests:**
            - Conversion between matrix and `neuroim2` objects preserves information
            - Spatial metadata is maintained throughout pipeline
            - Missing data (NaN) handling
        *   **Input validation tests:**
            - Mismatched dimensions
            - Invalid parameters (negative lambda, rank > min(V,K))
            - Non-finite values in Y
    *   Acceptance Criteria: Key risks mitigated in implementation. Robust handling of neuroimaging data characteristics and potential issues.

---

**🔄 SPRINT 1 INTEGRATION NOTES:**

**Shared Components Ready for CBD Reuse:**
- Project setup, dependencies, GitHub Actions (T01)
- Data structure utilities: matrix ↔ neuroimaging conversion (T02)  
- HRF utilities: `fmrireg` integration, convolution functions (T03)
- Simulation framework: unified data generation supporting both algorithms (T04)
- Visualization framework: spatial maps, convergence, state timecourses (T07)
- Testing infrastructure: fixtures, utilities, benchmarks (T14)
- Documentation patterns and package structure (T15)

**CLD → CBD Algorithm Differences:**
- **CLD:** GLM+SVD spatial learning, FISTA optimization, TV penalty
- **CBD:** VB parameter updates, HMM forward-backward, ELBO computation  
- **Shared:** Low-rank spatial patterns, HRF convolution, neuroimaging integration

**Sprint 1 Development Strategy:**
1. **Leverage shared infrastructure** from Sprint 1a without modification
2. **Implement CBD-specific algorithms** (VB updates, HMM) in parallel to CLD FISTA
3. **Extend simulation framework** to include Markov chain state generation
4. **Reuse visualization framework** with CBD-specific adaptations (ELBO plots, posterior uncertainty)
5. **Direct algorithm comparison** using shared simulation and benchmarking infrastructure

**Estimated Sprint 1 Efficiency Gain: ~40% reduction in development time through infrastructure reuse**

**Sprint 1a Summary:** This sprint delivers a working Continuous Linear Decoder prototype that integrates seamlessly with R's neuroimaging ecosystem. The implementation balances performance (via strategic use of RcppArmadillo and algorithmic optimizations) with usability (R6 interface, neuroimaging data structure support). The foundation is set for future enhancements including online learning, multi-subject analysis, and the full Bayesian treatment in subsequent sprints.

