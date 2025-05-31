Okay, building on the foundation laid in Sprint 1, Sprint 2 will focus on:
1.  **Introducing HRF Estimation:** Moving from fixed HRFs to estimating voxel-specific HRFs using a basis set.
2.  **Rcpp Optimizations:** Profiling Sprint 1 code and migrating identified bottlenecks (especially in the E-step likelihoods and forward-backward) to Rcpp.
3.  **Improving M-Step Updates:** Refining updates for `U`, `V`, and `H_v` (HRF coefficients).
4.  **Enhanced Diagnostics & Testing:** More rigorous testing and diagnostic plots.

**Sprint 2 Goal:** Implement voxel-specific HRF estimation within the VB framework. Significantly improve computational performance of critical E-step components by migrating them to Rcpp. Refine M-step updates and enhance model diagnostics. Continue building on established R neuroimaging infrastructure (`neuroim2`, `fmrireg`) while optimizing performance. **All new code must compile on macOS, Linux, and Windows (CRAN toolchains) with `R CMD check --as-cran` ≤ 1 NOTE.**

**Assumptions for Sprint 2:**
*   Sprint 1 deliverables (R6 class, basic VB loop, simulator) are functional and pass full package checks.
*   No GMRF prior on HRFs *yet* (spatial smoothness for HRFs will be Sprint 3). We'll estimate independent voxel HRFs first.
*   Use shared noise variance `σ²` initially – skip voxel-specific `σ_v²` until HRFs converge to avoid changing two major variance sources simultaneously.
*   Continue leveraging `fmrireg` HRF basis functions and `neuroim2` data structures.
*   **Continuous integration**: GitHub Actions matrix builds (macOS, ubuntu, windows) must stay green after each ticket merge.

---

**Sprint 2 Tickets:**

**Theme 1: HRF Estimation Implementation**

1.  **Ticket S2-T01: Define HRF Basis Set Functionality**
    *   Description: Create helper functions to generate HRF basis sets, building on `fmrireg`'s established HRF components (`HRF_SPMG1`, `HRF_SPMG2`, `HRF_SPMG3`, `HRF_FIR`, `HRF_BSPLINE`). The function should return a `L_h × L_basis` matrix.
        *   Integrate with `fmrireg`'s `hrf()` and basis construction functions where possible.
        *   Extend support for canonical HRF + temporal & dispersion derivatives; FIR basis sticks; or splines using `splines::bs`.
        *   Ensure compatibility with `fmrireg`'s sampling frame and time‐resolution conventions.
        *   **Pin TR & duration defaults** (e.g., TR = 0.8 s, duration = 32 s) so unit tests never drift.
        *   **Expose `orthonormal = TRUE` flag** to make later GMRF algebra nicer.
    *   Acceptance Criteria: Functions `create_hrf_basis_canonical(TR, duration_secs, n_derivatives)` and/or `create_hrf_basis_fir(TR, duration_secs, fir_resolution_secs)` are available and produce valid basis matrices compatible with `fmrireg` conventions. Unit tests cover TR drift edge cases.
2.  **Ticket S2-T02: Modify R6 `initialize` for HRF Parameters**
    *   Description: Update `ContinuousBayesianDecoder$initialize` to:
        *   Accept HRF basis type/parameters (compatible with `fmrireg` HRF specifications).
        *   **Initialize voxel‐specific HRF coefficients `private$.H_v` with least‐squares fit to each voxel's average response** (not zeros) – cuts EM iterations by 30-40 %.
        *   Maintain compatibility with both standard matrices and `neuroim2` data structures.
    *   Acceptance Criteria: Decoder initializes HRF coefficients optimally and stores the HRF basis. Initialization runtime < 30 s for 1 k voxels × 200 TR dataset.
3.  **Ticket S2-T03: Modify Simulator to Generate Data with Voxel-Specific HRFs**
    *   Description: Update `simulate_cbd_data` to accept `true_H_v` (coefficients) and an HRF basis (using `fmrireg` HRF functions for realistic HRF shapes). It should then generate `Y` using these voxel-specific HRFs.
        *   Each voxel `v` will have `h_v = B_HRF %*% true_H_v[v,]`.
        *   Convolution `(h_v ★ S_k,:)` now uses the voxel's specific `h_v`.
        *   Ensure compatibility with both standard matrix and `neuroim2::NeuroVec` output formats.
        *   **Return a structured list** with `NeuroVec`, `true_H`, `true_sigma2` to streamline unit tests.
    *   Acceptance Criteria: Simulator can generate data with known, varying HRFs across voxels using realistic `fmrireg` HRF basis functions. Works with both data structure types.

**Theme 2: Rcpp Optimization of E-Step (Before HRF M-Step)**

4.  **Ticket S2-T04: Profile Sprint 1 E-Step Code**
    *   Description: Use `profvis`, `lineprof`, and `Rprofmem` to identify performance & memory hotspots in `private$.compute_log_likelihoods`, `private$.forward_pass`, and `private$.backward_pass`.
    *   **Tag lines** with `# PERF:` comments to make future diffs obvious.
    *   Acceptance Criteria: Profiling markdown report committed in `inst/benchmarks`. Key bottleneck list produced.
5.  **Ticket S2-T05: Rcpp Implementation – Voxel-Specific HRF Convolution**
    *   **Minor tweak**: Expose `fft_threshold` arg so unit tests can force direct-time convolution.
6.  **Ticket S2-T06: Rcpp Implementation: Likelihood Calculation Kernel**
    *   Description: Migrate the core loop of `private$.compute_log_likelihoods` to Rcpp. This function will take `Y_block` (`V x T_block`), `U`, `V`, the pre-computed voxel-specific convolved state regressors (from S2-T05), and `sigma2` (shared for now).
        *   **Project Y once**: `Y_r = UᵀY`; keep both `Y_r` and `W_r = Σ_r Vᵀ` in Rcpp so GEMM is `r×T` not `V×T`. Big win if `V > 20k`.
        *   Handle data efficiently regardless of whether it originates from standard matrices or `neuroim2` structures.
    *   Acceptance Criteria: Rcpp function `compute_log_likelihoods_rcpp` significantly outperforms the R version and works with both data structure types.
7.  **Ticket S2-T07: Rcpp Implementation: Forward Algorithm**
    *   Description: Convert `private$.forward_pass` to Rcpp. Input: `log_likelihoods_rcpp` (from S2-T06), `Pi`, `pi0`. Output: `alpha_rcpp` matrix and data log-likelihood. Ensure log-sum-exp for numerical stability.
        *   **Store α/β in column-major Armadillo mats** and reuse in ξ calculation to cut copies.
    *   Acceptance Criteria: `forward_pass_rcpp` is faster and numerically stable.
8.  **Ticket S2-T08: Rcpp Backward Algorithm**
    *   **🔄** Ensure shared code path with forward implementation to avoid duplication (template or header reuse).
9.  **Ticket S2-T09: Integrate Rcpp Components & Engine Switch**
    *   **Add unit test** that toggles `engine = "R"/"cpp"` and asserts identical α/β within 1e-6.

**Theme 3: HRF M-Step Implementation (After Rcpp E-Step)**

10. **Ticket S2-T10: Batched HRF Coefficient Update**
    *   **Watch-out:** Armadillo `solve_sympd` can silently fallback to dense if the system is not SPD. Check `info` flag and fall back to `solve()` with regularisation if needed.

**Theme 4: M-Step Refinements**

11. **Ticket S2-T11: Modify Likelihood Calculation for Voxel-Specific HRFs**
    *   Description: Update `private$.compute_log_likelihoods` to use current estimates of `private$.H_v` and `private$.hrf_basis` to compute voxel-specific convolved regressors for each state when calculating `log p(Y_:,t | S_:,t = j, Θ_current)`.
        *   **Cache `H★S` per HRF basis column** before looping over states - you re-use those `T×L` slices many times in the likelihood kernel.
        *   Maintain compatibility with both standard matrices and `neuroim2` data structures.
    *   Acceptance Criteria: Log-likelihood calculation correctly incorporates current estimates of voxel-specific HRFs and works with both data structure types.

**Theme 5: Diagnostics, Testing & Documentation**

12. **Ticket S2-T12: Refine M-Step for `U` and `V` with Estimated HRFs**
    *   Description: Update the M-step for `U` and `V` to correctly incorporate `E_q[H_v]` (or current point estimates of `H_v`) and `E_q[S]`.
        *   **Do one SVD of W each iteration** then update `U`, `V` via Procrustes; otherwise alternating LS can diverge once HRFs start moving.
        *   Consider projecting `Y` onto `U`'s current estimate (`Y_proj = UᵀY`) and performing `V` updates in the low-rank space.
        *   Ensure efficient handling of both standard matrices and `neuroim2` data structures.
    *   Acceptance Criteria: Updates for `U` and `V` are robust and theoretically sound with estimated HRFs, working with both data structure types.
13. **Ticket S2-T13: Update ELBO Calculation for Estimated HRFs**
    *   Description: Modify the ELBO calculation to include terms for:
        *   `E_q[log p(Y | S, U, V, {H_v}, σ²)]` using voxel-specific HRFs and shared noise.
        *   Priors on `H_v` (e.g., simple Gaussian `N(0, λ_I⁻¹I)` for coefficients initially).
        *   `E_q[log p(H_v)] - E_q[log q(H_v)]`.
        *   **Keep "prior-minus-entropy" terms optional** behind `self$config$elbo_full`.
    *   Acceptance Criteria: ELBO calculation is comprehensive for the current model complexity.
14. **Ticket S2-T14: Diagnostic Plot for Estimated HRFs**
    *   **Minor tweak**: Add `geom_ribbon()` for 95 % CI if posterior variance available.
15. **Ticket S2-T15: Unit Tests for Rcpp E-Step Components**
    *   **Add Windows toolchain CI** (Rtools42) to ensure templates compile without `-fopenmp` if not available.
16. **Ticket S2-T16: End-to-End Test with Simulated Data (Voxel-Specific HRFs)**
    *   Description: Update the end-to-end test script:
        1.  Simulate data with known, varying `H_v` (using `fmrireg` HRF functions).
        2.  Fit the CBD model.
        3.  Check recovery of `H_v` (e.g., correlation between true and estimated `H_v` coefficients for some voxels).
        4.  Check ELBO convergence.
        5.  **Add timing benchmark**: `system.time()` before/after Rcpp swap; target ≥ 5× speed improvement.
        6.  Test with both standard matrices and `neuroim2` data structures.
    *   Acceptance Criteria: Model can recover HRFs to some degree. Noticeable performance improvement from Rcpp is observed (≥ 5× speedup). Both data structure types work correctly.
17. **Ticket S2-T17: Document Rcpp Functions and HRF Implementation Choices**
    *   Description: Add comments and documentation for new Rcpp functions. Document how HRF basis sets are handled and how `H_v` is estimated. Document integration with `fmrireg` HRF components and `neuroim2` data structures.
        *   **Note**: Any Rcpp function exporting Armadillo sparse mats must include `// [[Rcpp::depends(RcppArmadillo)]]` and `// [[Rcpp::plugins(cpp17)]]` to silence CRAN checks.
    *   Acceptance Criteria: New code is adequately documented, including integration with neuroimaging packages.
18. **Ticket S2-T18: Code Review & Refactor Sprint 2 Components**
    *   Description: Review all new and modified code from Sprint 2. Identify areas for further optimization or simplification. Plan for GMRF prior introduction in Sprint 3. Ensure proper integration with `neuroim2` and `fmrireg` components throughout.
    *   Acceptance Criteria: Sprint 2 code reviewed. Plan for next sprint refined. Integration with neuroimaging ecosystem verified.

---

**Key Ordering & Scope Notes:**
1. **Rcpp kernels (T04-T09) come before HRF M-step (T10-T11)** - need profiling numbers on the old M-step to know if HRF update is suddenly dominating.
2. **Skip voxel-wise `σ_v²` this sprint** - changing two major variance sources simultaneously can hide bugs. Accept shared `σ²` until unit tests on `H_v` pass.
3. **Memory management**: Convolved design tensor stored as `(K × T) × L_basis`, not `V × K × T`.

This sprint significantly increases model complexity (HRF estimation) and technical complexity (Rcpp). With these optimizations, you'll achieve credible HRFs, ~5× wall-time gain, and a clean runway for GMRF spatial priors in Sprint 3. It's crucial to test incrementally and ensure that the Rcpp components are both correct and faster while maintaining seamless integration with the established R neuroimaging infrastructure.

---

**Cross-Cutting Risks (Sprint 2):**
*   **FFT on Windows:** `stats::fft` call inside Rcpp uses R's FFTPACK (no OpenMP). Document possible slowdown.
*   **Sparse → dense spill-over:** Ensure big Laplacian matrices stay sparse during HRF update.
*   **Compiler flags:** Avoid `-march=native`; CRAN disallows.
*   **Memory duplication inside `as.matrix(NeuroVec)`:** Large ROI runs may duplicate data; consider `as(view(Y))` thin wrapper instead.

---

**Exit Criteria for Sprint 2:**
1.  Rcpp E-step passes equivalence tests (α/β max abs diff < 1e-6).
2.  ≥ 5× speed-up on 1 k vox × 200 TR synthetic dataset compared to pure R.
3.  HRF coefficients recover synthetic truth (median voxel corr > 0.7).
4.  GitHub Actions CI green on all OS.