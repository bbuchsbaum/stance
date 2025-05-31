Okay, let's break down Sprint 1 into granular tickets. The goal for Sprint 1 is to lay the foundational R6 class structure, implement the core mathematical components for a *simplified* version of the model (perhaps fixed HRF, no GMRF yet), and get the basic VB loop running with simulated data. We'll focus on correctness and structure before heavy Rcpp optimization, but keep Rcpp in mind for design.

**Sprint 1 Goal:** Establish a functional R6 class for the Continuous Bayesian Decoder (CBD) with a basic Variational Bayes (VB) inference loop, capable of running on simulated data with a simplified model (e.g., known/fixed HRFs, no GMRF, shared noise variance). Implement core HMM forward-backward in R. Build on established R neuroimaging infrastructure using `neuroim2` for data structures and `fmrireg` for fMRI regression components. **This sprint leverages ~60-70% of the shared infrastructure established in Sprint 1a (CLD), including project setup, neuroimaging data integration, HRF utilities, simulation framework, and testing infrastructure. Focus is on CBD-specific algorithm implementation while reusing the proven foundation.** **All code will follow standard R package conventions (roxygen2 doc, testthat, pkgdown site).**

**Assumptions for Sprint 1:**
*   **Sprint 1a foundation completed:** All shared infrastructure from CLD sprint is available and tested
*   Focus on a single, known HRF initially, or a very simple HRF basis (e.g., just canonical from `fmrireg`).
*   Shared isotropic noise variance (`σ²`) across all voxels.
*   No GMRF prior on HRFs yet (as HRFs are fixed or very simple).
*   Low-rank factorization `W = UVᵀ` is a core component from the start (where U is V×r, V is K×r, and W is V×K of rank r).
*   Focus on R implementation first, identifying Rcpp candidates for Sprint 2.
*   Leverage shared utilities: `neuroim2` integration, `fmrireg` HRF components, simulation framework.
*   **Project should R CMD check cleanly (≤1 NOTE) on all OSes.**

---

**Sprint 1 Tickets:**

**Theme 1: Leveraging Sprint 1a Infrastructure (🟢 REUSED COMPONENTS)**

1.  **Ticket S1-T01: Validate & Extend Shared Infrastructure** 🟢 **[BUILDS ON SPRINT 1a]**
    *   Description: Verify Sprint 1a shared infrastructure works for CBD requirements:
        *   Test `R/data_structures.R` utilities with CBD-specific data flows
        *   Validate `R/hrf_utils.R` supports HMM convolution operations needed for CBD
        *   Extend `R/simulate.R` to include Markov chain state generation for CBD validation
        *   Check `R/visualization.R` framework can accommodate ELBO convergence plots
        *   **🔄 No major changes expected** - mainly validation and minor extensions
    *   Acceptance Criteria: Sprint 1a infrastructure validated for CBD use. Minor extensions completed. **Foundation ready for CBD algorithm implementation.**
2.  **Ticket S1-T02: Extend Simulation Framework for CBD** 🟢 **[BUILDS ON SPRINT 1a-T04]**
    *   Description: Extend `simulate_fmri_data()` from Sprint 1a to fully support CBD mode:
        *   Add Markov chain state sequence generation with transition matrix Π
        *   Support low-rank spatial pattern generation `W = UVᵀ`
        *   Include HMM parameter ground truth for CBD validation
        *   Maintain compatibility with CLD mode for algorithm comparison
        *   **🔄 Reuse existing infrastructure** from Sprint 1a with CBD-specific additions
    *   Acceptance Criteria: Unified simulation function supports both CLD and CBD. **Enables direct algorithm comparison and CBD parameter recovery validation.**

**Theme 2: CBD-Specific R6 Class Structure**

3.  **Ticket S1-T03: Define `ContinuousBayesianDecoder` R6 Class Skeleton** 🟡 **[NEW - FOLLOWS SPRINT 1a PATTERNS]**
    *   Description: Create `R/continuous_bayesian_decoder.R` following Sprint 1a R6 patterns:
        *   Private fields: `$.Y_data`, `$.Y_proj` (UᵀY), `$.U`, `$.V`, `$.Pi`, `$.pi0`, `$.sigma2`, `$.S_gamma`, `$.S_xi`, `$.hrf_kernel`, `$.elbo_history`, `$.config` (list with max_iter, tol, rank, n_states), `$.neuro_metadata` (for NeuroSpace if applicable)
        *   Public methods: `initialize()`, `fit()`, `get_spatial_maps()`, `get_state_posteriors()`, `get_parameters()`, `predict()`, `plot_convergence()`
        *   Private methods: `$.compute_log_likelihoods()`, `$.forward_backward()`, `$.update_U_V()`, `$.update_Pi()`, `$.update_sigma2()`, `$.compute_elbo()`
        *   **🔄 Reuse Sprint 1a data structure patterns** for `neuroim2` integration
    *   Acceptance Criteria: R6 class definition is syntactically correct. roxygen2 skeleton added. Handles both matrix and `neuroim2::NeuroVec` input types. **Follows proven Sprint 1a design patterns.**
4.  **Ticket S1-T04: Implement Robust `initialize` Method** 🟡 **[NEW - USES SPRINT 1a UTILITIES]**
    *   Description: 
        *   **Input handling**: Use `validate_fmri_input()` and `extract_data_matrix()` from Sprint 1a
        *   **Validation**: Check dimensions, non-finite values, ensure T > K using shared validation utilities
        *   **HRF setup**: Use `setup_hrf_kernel()` from Sprint 1a HRF utilities
        *   **Initialization**:
            - `U`: Initialize via truncated SVD of Y: `svd(Y, nu = rank, nv = 0)$u`
            - `V`: Random normal, then orthogonalize columns via QR
            - `Pi`: Slightly diagonal-dominant: `diag(K) * 0.8 + matrix(0.2/K, K, K)`
            - `pi0`: Uniform `rep(1/K, K)`
            - `sigma2`: `0.1 * var(as.vector(Y))`
        *   **Precompute**: `Y_proj = crossprod(U, Y)` for efficiency
        *   **🔄 Leverage all Sprint 1a infrastructure** for data handling and validation
    *   Acceptance Criteria: Decoder instantiates in <1s for V=10000, T=500. All parameters have correct dimensions. Passes validation tests. **Uses shared utilities seamlessly.**
5.  **Ticket S1-T05: Extend Shared S3 Methods for CBD** 🟢 **[BUILDS ON SPRINT 1a-T07]**
    *   Description: 
        *   Extend `R/visualization.R` from Sprint 1a with CBD-specific methods:
            - `print.ContinuousBayesianDecoder()`: Show dimensions, convergence status, final ELBO
            - `summary.ContinuousBayesianDecoder()`: Detailed parameter summary with posterior state probabilities  
            - Extend `plot_convergence()` to handle ELBO traces (vs FISTA objective)
            - `coef.ContinuousBayesianDecoder()`: Extract key parameters as named list
            - `fitted.ContinuousBayesianDecoder()`: Return fitted values Ŷ = W(H★S_gamma)
        *   **🔄 Reuse Sprint 1a visualization framework** with CBD-specific adaptations
    *   Acceptance Criteria: All S3 methods documented and tested. CBD plots work alongside CLD plots from Sprint 1a. **Unified visualization framework for algorithm comparison.**

**Theme 3: VB Loop – Pure R Prototype**

6.  **Ticket S1-T06: Implement Likelihood Calculation** 🟡 **[NEW ALGORITHM - USES SPRINT 1a HRF UTILITIES]**
    *   Description: Create `private$.compute_log_likelihoods()`:
        *   Work in projected space: Use pre-computed `Y_proj = UᵀY` (r×T)
        *   For each state k and time t, compute predicted signal in low-rank space:
            - Create state indicator for state k active at time t
            - Convolve with HRF using `convolve_with_hrf()` from Sprint 1a HRF utilities
            - Project through V: `mu_proj_k = V[k,] * (h ★ s_k)`
        *   Log-likelihood: `-0.5/sigma2 * ||Y_proj[,t] - mu_proj_k||^2 - r*log(2*pi*sigma2)/2`
        *   Use log-sum-exp trick for numerical stability when needed
        *   **🔄 Leverage Sprint 1a HRF convolution infrastructure**
    *   Acceptance Criteria: Returns K×T matrix. Likelihoods increase when true state is active. Numerically stable for sigma2 → 0. **Uses shared HRF utilities.**
7.  **Ticket S1-T07: HMM Forward-Backward Algorithm** 🟡 **[NEW ALGORITHM]**
    *   Description: Implement `private$.forward_backward(log_likelihoods, Pi, pi0)`:
        *   **Forward pass with scaling**:
            - α[1,] = pi0 * exp(log_lik[,1]); c[1] = 1/sum(α[1,])
            - For t=2:T: α[t,] = (α[t-1,] %*% Pi) * exp(log_lik[,t]) * c[t]
            - Track log-likelihood: sum(log(1/c))
        *   **Backward pass**:
            - β[T,] = 1
            - For t=(T-1):1: β[t,] = Pi %*% (β[t+1,] * exp(log_lik[,t+1])) * c[t+1]
        *   **Posteriors**:
            - γ[t,k] = α[t,k] * β[t,k] (already normalized)
            - ξ[t,i,j] ∝ α[t,i] * Pi[i,j] * exp(log_lik[j,t+1]) * β[t+1,j]
        *   Store γ in `private$.S_gamma`, ξ in `private$.S_xi`
    *   Acceptance Criteria: All posteriors sum to 1. No NaN/Inf values. Matches HMM reference implementation. **Pure CBD algorithm - no shared components.**
8.  **Ticket S1-T08: VB M-Step Parameter Updates** 🟡 **[NEW ALGORITHM]**
    *   Description: Implement parameter update methods:
        *   **`private$.update_U_V()`**: 
            - Define expected signal: X_expected = H ★ S_gamma (using current posteriors)
            - Update V: For each k, `V[k,] = solve(crossprod(U) + lambda*I) %*% crossprod(U, Y %*% t(X_expected[k,]))`
            - Update U: `U_new = Y %*% t(X_expected) %*% solve(crossprod(X_expected) + lambda*I)`
            - Orthonormalize U via QR: `U = qr.Q(qr(U_new))`
            - Recompute Y_proj = crossprod(U, Y)
        *   **`private$.update_Pi()`**: 
            - `Pi[i,j] = (sum_t ξ[t,i,j] + prior_pseudocount) / (sum_t γ[t,i] + K*prior_pseudocount)`
            - Ensure rows sum to 1
        *   **`private$.update_sigma2()`**:
            - Compute expected residuals using current parameters
            - `sigma2 = mean((Y - W*(H★S_gamma))^2)`
    *   Acceptance Criteria: Parameters converge toward true values on simulated data. U remains orthonormal. **CBD-specific algorithms.**
9.  **Ticket S1-T09: ELBO Calculation & Convergence** 🟡 **[NEW ALGORITHM]**
    *   Description: Implement `private$.compute_elbo()`:
        *   **Expected log-likelihood**: Use current log_likelihoods and γ
        *   **HMM entropy**: `-sum(γ * log(γ + 1e-10))`
        *   **Prior terms**: `sum(log(pi0) * γ[1,])` + transition terms from ξ and Pi
        *   **Parameter priors** (if using proper Bayesian): Gaussian for U,V; Dirichlet for Pi
        *   Track history in `private$.elbo_history`
        *   Convergence: `abs(elbo_new - elbo_old) / abs(elbo_old) < tol`
        *   **🔄 Extend Sprint 1a convergence plotting** to handle ELBO traces
    *   Acceptance Criteria: ELBO increases monotonically (allowing small numerical errors ~1e-6). Converges on test data. **Uses shared visualization infrastructure.**

**Theme 4: Integration & Validation**

10. **Ticket S1-T10: Main `fit()` Method** 🟡 **[NEW - FOLLOWS SPRINT 1a PATTERNS]**
    *   Description: Implement public `fit(max_iter = 100, tol = 1e-4, verbose = TRUE, save_history = FALSE)` following Sprint 1a patterns:
        *   Initialize progress bar with `cli::cli_progress_bar(total = max_iter)`
        *   Main VB loop:
            ```r
            for (iter in 1:max_iter) {
              # E-step
              log_lik <- private$.compute_log_likelihoods()
              private$.forward_backward(log_lik, self$Pi, self$pi0)
              
              # M-step  
              private$.update_Pi()
              private$.update_U_V()
              private$.update_sigma2()
              
              # Convergence
              elbo <- private$.compute_elbo()
              if (private$.check_convergence(elbo, tol)) break
              
              cli::cli_progress_update()
            }
            ```
        *   Store convergence status, final ELBO, iteration count
        *   Return self invisibly for method chaining
        *   **🔄 Follow Sprint 1a interface patterns** for consistency
    *   Acceptance Criteria: Fits V=1000, T=200 model in <20s. Progress bar shows ETA. Convergence detected properly. **Consistent with Sprint 1a CLD interface.**
11. **Ticket S1-T11: CBD-Specific Testing** 🟢 **[BUILDS ON SPRINT 1a-T14]**
    *   Description: 
        *   Extend `tests/testthat/test-infrastructure.R` with CBD-specific tests
        *   Create `tests/testthat/test-cbd.R` for CBD algorithm validation:
            - HMM forward-backward accuracy vs reference implementations  
            - VB parameter update convergence on simulated data
            - ELBO monotonicity and convergence
            - Parameter recovery metrics using shared simulation framework
        *   **🔄 Reuse Sprint 1a test fixtures and utilities**
        *   Test CBD vs CLD on identical simulated datasets for algorithm comparison
    *   Acceptance Criteria: Comprehensive CBD test suite. **Builds on Sprint 1a testing infrastructure.** Algorithm comparison tests available.
12. **Ticket S1-T12: Performance Profiling & Sprint 2 Planning** 🟡 **[NEW - FOLLOWS SPRINT 1a PATTERNS]**
    *   Description: 
        *   Profile full VB fit on V=5000, T=500 dataset using `profvis`
        *   Compare CBD vs CLD performance using shared benchmarking infrastructure from Sprint 1a
        *   Identify bottlenecks (expected: convolution, forward-backward, matrix multiplies)
        *   Create `inst/profiling/cbd_profile_results.md` documenting:
            - Time breakdown by function
            - Memory allocation patterns  
            - Recommendations for Rcpp conversion priority
            - Performance comparison with Sprint 1a CLD
        *   Draft Sprint 2 plan focusing on CBD-specific optimizations
    *   Acceptance Criteria: Profiling identifies >80% of runtime. **Algorithm comparison available.** Sprint 2 plan includes specific Rcpp targets for CBD.

**Theme 5: Documentation & Integration**

13. **Ticket S1-T13: Extend Shared Documentation** 🟢 **[BUILDS ON SPRINT 1a-T15]**
    *   Description: 
        *   Extend `vignettes/getting-started.Rmd` to include CBD examples alongside CLD
        *   Update `README.md` with algorithm comparison and selection guidance
        *   Add CBD examples to `inst/examples/` directory
        *   Document when to use CBD vs CLD based on use case requirements
        *   **🔄 Leverage Sprint 1a documentation infrastructure**
        *   Update package-level documentation to describe both algorithms
    *   Acceptance Criteria: Comprehensive documentation covering both algorithms. **Clear algorithm selection guidance.** Builds on Sprint 1a documentation foundation.

---

**🔄 SPRINT 1A DEPENDENCY SUMMARY:**

**🟢 Directly Reused Components (No Changes):**
- Project setup, dependencies, GitHub Actions (S1a-T01)
- Data structure utilities: `validate_fmri_input()`, `extract_data_matrix()`, etc. (S1a-T02)  
- HRF utilities: `setup_hrf_kernel()`, `convolve_with_hrf()`, etc. (S1a-T03)
- Testing infrastructure: fixtures, utilities, benchmarks (S1a-T14)
- Package structure and dependency management (S1a-T15)

**🟡 Extended Components (Minor Changes):**
- Simulation framework: Add Markov chain generation to existing `simulate_fmri_data()` (S1a-T04 → S1-T02)
- Visualization framework: Add ELBO plotting to existing convergence plots (S1a-T07 → S1-T05)
- Documentation: Extend existing vignettes and README (S1a-T15 → S1-T13)

**🔴 New CBD-Specific Components:**
- R6 class for CBD: New algorithm implementation (S1-T03)
- VB algorithms: HMM forward-backward, parameter updates, ELBO (S1-T06 through S1-T09)
- Main fit method: VB loop implementation (S1-T10)

**Development Efficiency:**
- **~60-70% infrastructure reuse** reduces Sprint 1 development time
- **Shared simulation enables direct algorithm comparison** 
- **Common testing and benchmarking infrastructure** ensures consistent validation
- **Unified documentation** provides users with clear algorithm selection guidance

**🟢 PARALLEL EXECUTION STRATEGY (Updated):**

**Phase 1 (Validation):** S1-T01 → S1-T02 (validate and extend Sprint 1a infrastructure)

**Phase 2 (Parallel Development):**
- **🟢 Group A (Interface & Testing):** S1-T03, S1-T04, S1-T05, S1-T11, S1-T13 - Build on Sprint 1a patterns
- **🟢 Group B (VB Algorithms):** S1-T06, S1-T07, S1-T08, S1-T09 - Pure CBD algorithms  

**Phase 3 (Integration):** S1-T10 → S1-T12

**Merge Conflict Risk: SIGNIFICANTLY REDUCED** due to Sprint 1a foundation providing stable interfaces and shared utilities.

---

**Exit Criteria for Sprint 1:**
1.  All tickets S1-T01 through S1-T13 completed and merged to main
2.  `R CMD check` passes with 0 errors, 0 warnings, ≤1 NOTE on all platforms (via GitHub Actions)
3.  Test coverage ≥80% (via covr/codecov)
4.  Successful parameter recovery on simulated data:
    - Spatial correlation(U_true, U_est) > 0.7
    - State detection AUC > 0.85
    - Convergence in <50 iterations for standard problems
5.  Documentation complete:
    - All exported functions have roxygen2 docs with examples
    - Vignette builds and is informative
    - pkgdown site deploys successfully
6.  Performance benchmarks documented:
    - V=1000, T=200 fits in <20s
    - Memory usage <2GB for typical problems
    - Profiling report identifies optimization targets

**Sprint 1 Deliverables Summary:**
- Fully functional R package 'stance' implementing basic CBD
- Comprehensive test suite and documentation
- Clear roadmap for Sprint 2 optimizations

