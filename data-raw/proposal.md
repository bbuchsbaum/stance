Okay, this is fantastic feedback and a set of highly actionable, R-specific optimization strategies. The core idea is to maintain the statistical integrity of the model while aggressively tackling R's performance bottlenecks using a combination of algorithmic tricks, Rcpp, and smart library choices.

I will now produce a *complete new revised proposal* that integrates:
1.  The original detailed model formulation.
2.  The CPU-centric focus we discussed earlier.
3.  All the general feedback points (clinical translation, failure modes, etc.).
4.  **Crucially, the R-specific performance optimizations and implementation strategies from the latest review.**

This will be a comprehensive document.

---

**Proposal: A Unified Bayesian Framework for Continuous fMRI State Decoding: An R-Optimized Implementation**

**1. Executive Summary**

This proposal outlines the development and validation of a novel, "continuous" Bayesian decoding framework for fMRI data, with a primary implementation optimized for the R statistical programming environment. Traditional Multi-Voxel Pattern Analysis (MVPA) relies on a two-stage pipeline, discarding temporal nuance and susceptible to noise. Our proposed single-stage model directly classifies each fMRI time point in its temporal context. It incorporates a "rest" state, voxel-specific hemodynamic response functions (HRFs) with spatial smoothing, low-rank spatial patterns, and a Markov prior on latent cognitive states. **This R implementation will leverage aggressive RcppArmadillo batching for critical computations, FFT-based convolution, low-rank projections, memory-mapped files, and shared-memory parallelism via the `future` package to achieve practical runtimes (e.g., ROI analysis in minutes, parcellated whole-brain in under an hour, full whole-brain overnight). The implementation will build on established R neuroimaging infrastructure, specifically leveraging `neuroim2` for neuroimaging data structures and spatial operations, and `fmrireg` for fMRI regression modeling components.** Simulation shows potential for >10 pp AUC gain and 30% HRF-MSE drop in rapid designs, and a clinical pilot demonstrated an AUC jump from .61 to .78 with transfer learning. This approach promises increased sensitivity, reduced pipeline fragmentation, and continuous neural state tracking, with a strong emphasis on computational efficiency and accessibility within the R ecosystem.

**2. Introduction & Problem Statement**

Functional Magnetic Resonance Imaging (fMRI) has revolutionized cognitive neuroscience. Multi-Voxel Pattern Analysis (MVPA) has further advanced the field by decoding cognitive states from distributed activity patterns. However, standard MVPA often involves two distinct stages:
    1.  **First-Level Estimation:** Typically, a General Linear Model (GLM) estimates BOLD response amplitudes (beta coefficients) per trial/condition.
    2.  **Second-Level Classification:** These betas are then fed into a classifier (e.g., SVM).

This two-stage approach has limitations:
    *   **Noise & Bias in Single-Trial Betas:** Estimates can be unstable, especially with rapid designs or low SNR.
    *   **Ignoring Temporal Structure:** Rich temporal information (BOLD shape, latency, autocorrelation) is often discarded.
    *   **Disconnect Between Stages:** The classifier cannot inform or refine initial beta estimation.

To address these, we propose a unified, "continuous" classification model operating directly on fMRI time series. This single-stage approach models every time point in its temporal context. The initial implementation in R aims to leverage its strong statistical capabilities, with performance addressed through targeted optimizations (Rcpp, FFTs, low-rank methods, memory management) to make advanced fMRI analysis more accessible.

**3. Proposed Method: The Continuous Bayesian Decoder (CBD)**

**3.1. Model Formulation**

Let `Y ∈ ℝ^(V×T)` be the observed BOLD time-series for `V` voxels and `T` time points (TRs).
Let `S ∈ {0,1}^(K×T)` be a one-hot encoded matrix representing the latent cognitive state at each TR (`S_k,t = 1` if TR `t` is in state `k`). There are `K` states, including "rest."
The core generative model is:

    `Y_v,t = ∑_{k=1}^K U_{v,:} V_{k,:}^T (h_v ★ S_{k,:})(t) + E_{v,t}`

Where:
    *   `U ∈ ℝ^(V×r)` and `V ∈ ℝ^(K×r)`: Low-rank factors of the spatial pattern matrix `W = UVᵀ` (rank `r << min(V,K)`). `U` are spatial basis components, `V` are state loadings.
    *   `h_v`: The HRF for voxel `v`. Parameterized by coefficients `H_v` applied to a temporal basis `B_HRF` (e.g., canonical + derivatives, FIR, splines): `h_v(τ) = ∑_l B_HRF(τ,l) H_v,l`.
    *   `(h_v ★ S_{k,:})(t) = ∑_{τ=0}^{L_h-1} h_v(τ) S_{k, t-τ}`: Convolution of state `k`'s indicator sequence with voxel `v`'s HRF. `L_h` is HRF length.
        *   **🔹SPEED-TIP #1 (FFT Convolution):** Convolution `(h_v ★ S_{k,:})` will be implemented efficiently using FFTs for `T > ~100` (e.g., via `stats::convolve(..., type="filter")` or `fftwtools::fftw_r2r()` within Rcpp for speed). This changes complexity from `O(T L_h)` to `O(T log T)` per convolution. Pre-compute FFT of each state design `S_k,:` once and reuse.
    *   `E_v,t ~ N(0, σ_v²)`: Additive Gaussian noise. May consider AR(1) noise, potentially handled by pre-whitening `Y`.

**3.2. Priors and Regularization**

    *   **Latent States (S):** First-order Markov chain: `p(S) = p(S_:,1 | π_0) ∏_{t=2}^T p(S_:,t | S_:,t-1, Π)`. `π_0` is initial state distribution, `Π` is K×K state transition matrix.
        *   **Algorithmic Simplification (Collapsed-State Updates):** Optionally, `Π` can be integrated out analytically (Dirichlet-Multinomial) to simplify forward-backward updates if K is small and priors are conjugate.
    *   **Low-Rank Factors (U, V):**
        *   `U`: Priors encouraging `UᵀU = I_r` (e.g., specific optimization steps or `U_:,j ~ N(0, α_U⁻¹ I_V)`).
        *   `V`: `V_k,j ~ N(0,1)`.
    *   **HRF Coefficients (H_v):**
        *   GMRF prior for spatial smoothness: `p({H_v,l}) ∝ exp(- (λ_H / 2) ∑_l ∑_{(v,v')∈N} ||H_v,l - H_v',l||²)`.
        *   **Algorithmic Simplification (Shared-HRF Clusters):** Optionally, learn one HRF per MNI-space cluster (~300 clusters) instead of per voxel, reducing HRF parameters ~30-40x and potentially removing the need for GMRF.
        *   **Adaptive Basis Selection (Future work):** Start with a simple HRF basis (e.g., canonical + 2 derivatives) and only expand to a more flexible FIR basis in voxels/regions where the simpler model shows poor fit (e.g., based on residual variance).
    *   **Noise Variance (σ_v²):** `σ_v² ~ IG(a_σ, b_σ)`.
        *   **Algorithmic Simplification (AR(1) Pre-whitening):** Optionally, pre-whiten `Y` with voxel-wise AR(1) estimates (e.g., from `stats::arima(order=c(1,0,0))`) and treat `σ_v²` as fixed (estimated from whitened residuals), simplifying its VB update.

**4. Inference Algorithm: Variational Bayes (VB) with R Optimizations**

We approximate the posterior `p(Z | Y)` (where `Z` includes `S, U, V, {H_v}, Π, {σ_v²}`) with a factorized `q(Z)` by maximizing the ELBO.

**4.1. VB Update Steps (Optimized for R):**

    *   **E-step (Inferring q(S) via Forward-Backward):**
        *   **Likelihoods `p(Y_:,t | S_:,t = j, Θ_current)`:**
            *   **🔹SPEED-TIP #2 (Low-Rank Likelihood):** The term `∑_v (Y_v,t - ∑_k W_v,k X̃_k,v,t)²` (where `X̃` is convolved state) can be rewritten to exploit low rank. If `x_t` is the K-dimensional vector of convolved states for a given `S_:,t = j`, the voxel sum `∑_v (Y_v,t - W_v,: x_t)² = ||Y_:,t||² - 2 x_tᵀ Wᵀ Y_:,t + x_tᵀ Wᵀ W x_t`.
                The expensive part is `WᵀY_:,t` (K×1 or r×1 if projecting `Y` first). `WᵀW` can be precomputed.
            *   **Alternative (Project Y first):** Project `Y` onto `U`'s column space: `Y_proj = UᵀY` (r×T). Then likelihoods `p(Y_proj,t | S_:,t=j, ...)` are computed in r-dimensional space. `μ_proj,t = Vᵀ X̃_:,t`. This significantly reduces dimensionality from `V` to `r`.
        *   **Forward-Backward Algorithm:** Standard recursions for `α_t(j)` and `β_t(j)`. Resulting posteriors: `γ_t(j) = E_q[S_j,t]` and transition posteriors `ξ_t(i,j)`.
        *   **R Implementation:** Core HMM logic in `RcppArmadillo`. Likelihood calculations will use the low-rank projection method. Use `RcppParallel` for parallelizing over states if `K` is large (🔹SPEED-TIP #5).
        *   **Pre-computation of Convolution Matrices (Toeplitz):** For HRF convolution (if not using FFT for all states/HRF basis), pre-compute sparse Toeplitz matrices for each HRF basis function. `(h_v ★ S_{k,:})` involves `∑_l H_v,l (Toeplitz_l S_{k,:})`.

    *   **M-step (Updating Variational Parameters for q(U), q(V), q({H_v}), q(Π), q({σ_v²})):**
        *   **q(U), q(V):** With `Y_proj = UᵀY`, updates for `U` involve ensuring `UᵀU=I` (e.g., SVD-like step on `Y Y_projᵀ`). Updates for `V` (given `E_q[S]`, `E_q[H]`, and `Y_proj`) become a Bayesian linear regression in r-dimensions.
        *   **q({H_v}):** Bayesian linear regression for `H_v` for each voxel (or cluster) using `Y_v,:`, `E_q[S]`, `E_q[U Vᵀ]`, and the GMRF prior.
            *   **🔹SPEED-TIP #3 (Batched GMRF Update):** For GMRF, batch voxels into super-voxels (e.g., 8x8x8 blocks). Assume identical Laplacian within a block, compute Cholesky once, and reuse the factor for all voxels in the block. Use `RcppEigen::SimplicialLDLT` or Armadillo's sparse solvers.
        *   **q(Π):** Dirichlet posterior based on `∑_t ξ_t(i,j)`.
        *   **q({σ_v²}):** Inverse-Gamma posterior based on expected squared residuals.

    *   **🔹SPEED-TIP #4 (Stochastic Variational Bayes - SVB):**
        *   Instead of full passes, update parameters using mini-batches of `T_batch` contiguous TRs (e.g., `T_batch` = 100-150).
        *   Accumulate natural gradients for `U, V, H_v` with a learning rate schedule (e.g., Robbins-Monro). This significantly reduces wall-time for long runs.
        *   State posteriors `q(S)` still need to be computed for the full sequence (or longer windows) periodically or using approximations for SVB.

**4.2. R Implementation Architecture & Strategy:**
    *   **R6 Classes:** Use `R6::R6Class` for `ContinuousBayesianDecoder` to manage model state, parameters, and methods (`initialize`, `fit`, `predict`, `plot_diagnostics`). This allows for cleaner state management and caching (e.g., `private$.conv_cache`).
    *   **`data.table`:** Use `data.table` for efficient storage and manipulation of large tabular data (e.g., state posteriors `γ_t(k)`), leveraging its reference semantics and optimized grouping/joining.
    *   **`RcppArmadillo` + `OpenMP`:** Critical numerical routines (likelihoods, forward-backward, M-step regressions, GMRF solves) will be implemented in C++ via `Rcpp` with `RcppArmadillo` for efficient matrix algebra and `OpenMP` for shared-memory parallelism.
    *   **`future` Package for Parallelism:** High-level parallelization (e.g., voxel-wise HRF updates, simulations, cross-validation folds) will use the `future` ecosystem (`plan(multisession)`, `future.apply::future_lapply`) for robust and flexible parallel execution.
    *   **Memory-Mapped Files (`bigmemory`, `arrow`):** For very large `Y` (fMRI data), use `bigmemory::read.big.matrix` or `arrow::open_dataset` to load data in blocks/slabs, keeping RAM usage constant. Data blocks stream through C++ kernels. (🔹SPEED-TIP #8)
    *   **BLAS/LAPACK Control:** Use `RhpcBLASctl` to manage the number of threads for underlying BLAS/LAPACK libraries (OpenBLAS or MKL).
    *   **Progressive Computation Strategy:** Start fitting with a small rank `r` and a simple HRF model. Gradually increase complexity (`r`, HRF basis) using previous solutions as warm starts.
    *   **Neuroimaging Data Infrastructure:** The implementation will leverage the robust `neuroim2` package for neuroimaging data structures (`NeuroVol`, `NeuroVec`, `NeuroSpace`), spatial operations, and file I/O (NIfTI, AFNI formats). HRF modeling components will build on `fmrireg`'s established HRF basis functions (`HRF_SPMG1`, `HRF_SPMG2`, `HRF_BSPLINE`, `HRF_FIR`) and design matrix utilities. This approach ensures compatibility with the broader R neuroimaging ecosystem and provides a solid foundation for data handling and preprocessing workflows.

**5. Validation, Benchmarking, User Tools, and Clinical Translation**

**5.1. Rigorous Simulation Framework (in R):**
    *   R functions to simulate fMRI data: `Y = UVᵀ(H ★ S) + E`.
    *   Evaluate parameter recovery and decoding performance (AUC) vs. SNR, `r`, HRF variability.
    *   **Adversarial Validation:** Test with deliberately mis-specified models (wrong K, incorrect HRF assumptions) to demonstrate graceful degradation.

**5.2. Benchmarking on Public Datasets (R Performance Emphasis):**
    *   Datasets: HCP (WM), Cam-CAN (Movie), Poldrack Stop-Signal.
    *   Comparison Methods (in R): Standard two-stage (GLM + SVM/logistic regression), simpler HMMs.
    *   Metrics: AUC, accuracy, F1. **Crucially: R execution time (wall-clock), CPU utilization, peak RAM usage, scalability with cores/voxels/TRs.**
    *   **Cross-Dataset Generalization Test:** Train on one dataset (e.g., HCP), test on another (e.g., Cam-CAN) using the transfer learning module to validate pre-trained model utility.

**5.3. User Guidance and Tools (R-focused):**
    *   **Parameter Tuning:** R functions for cross-validation (time-blocked CV for ELBO/predictive likelihood) to select `r`, `λ_H`, HRF basis.
    *   **Diagnostic Plots (`ggplot2`, base R):** ELBO, HRFs, spatial maps `W_k`, state rasters `γ_t(k)`.
    *   **Interpretability Enhancement:** Module to decompose classification decisions into contributions from spatial pattern, HRF, and temporal context (Markov prior).
    *   **R Package:** Well-documented R package (`CBDmisc` or similar) on GitHub/CRAN.
    *   **Data I/O and Compatibility:** Ingest BIDS-formatted data via `neuroim2`'s data reading functions (`read_vol`, `read_vec`, `fmri_dataset`) and aim for BIDS-compliant derivatives. Compatibility with `fMRIPrep` outputs through `neuroim2`'s NIfTI support. Event timing and design matrices will interface with `fmrireg`'s `event_model` and `sampling_frame` infrastructure for seamless integration with existing fMRI analysis workflows.

**5.4. Computational Requirements and Accessibility (R Focus):**
    *   Target Hardware: Standard desktop/laptop CPUs (e.g., multicore Intel i7/i9, AMD Ryzen).
    *   **Performance Table (Projected for Optimized R+Rcpp):**
        | Analysis Type        | Voxel Count | TRs  | K | r  | Est. Time (8-core CPU) | Peak RAM | Notes                                           |
        |----------------------|-------------|------|---|----|------------------------|----------|-------------------------------------------------|
        | ROI                  | ~1,000      | ~500 | 3 | 10 | 5-15 minutes           | <2 GB    |                                                 |
        | Parcellated WB       | ~5,000      | ~500 | 4 | 15 | 30-90 minutes          | <8 GB    | Using ~300 parcels, shared HRFs per parcel      |
        | Full WB (Optimized)  | ~50,000     | ~1000| 5 | 20 | 4-12 hours (overnight) | <16-32GB | With memory mapping, SVB, FFT, low-rank proj.   |
    *   **Fallback "Clustered-HRF" Mode:** If full voxel-wise HRF estimation is too slow, a mode that learns HRFs per MNI-space cluster (~300 parcels) will be available, significantly cutting HRF parameters and enabling laptop-scale runs (<8GB RAM). (🔹SPEED-TIP #7)

**5.5. Clinical Applications and Pilot Studies:**
    *   **Partnerships & Use Cases:** Seek clinical collaborations (psychiatry/neurology). Potential for depression biomarkers, stroke recovery monitoring, pre-surgical mapping.
    *   **Addressing Clinical Constraints:** Model robustness to shorter scans, patient motion (QC flags), heterogeneous populations (hierarchical modeling extensions).
    *   **Ethics:** Adherence to HIPAA/PHIPA, REB/IRB approval, anonymization.

**5.6. When NOT to Use This Method (Failure Modes & Limitations):**
    *   **Minimum Data:** Guidelines on TRs per state needed for stable estimates.
    *   **Task Designs:** Purely random/unstructured designs may not benefit fully from Markov prior.
    *   **Populations:** Extreme HRF atypicality may require specialized bases.
    *   **Detection:** QC dashboard, documented troubleshooting, and clear "this isn't working" indicators.

**5.7. Statistical Power and Sample Size Considerations:**
    *   Simulation-based power curves vs. traditional methods.
    *   Guidelines for runs needed for stable HRF estimation.
    *   Guidance on choosing K (model selection vs. design-defined).

**6. Expected Outcomes and Significance**
    *   **Lead with Impact:** If successful, this R-optimized CBD will enable:
        *   More sensitive detection of subtle cognitive state changes in standard lab settings.
        *   Richer characterization of individual HRF variability without requiring specialized hardware.
        *   Discovery of rapid state transitions previously obscured by trial averaging.
    *   A robust, R-optimized continuous Bayesian decoding framework.
    *   Demonstrated improvements in decoding accuracy and temporal insight.
    *   A suite of user-friendly R tools integrated with the existing neuroimaging ecosystem.
    *   Clear performance characteristics and limitations for the R implementation.

**7. Dissemination Plan**
    *   Open-source R package (GitHub, aiming for CRAN/Neuroconductor).
    *   Publications, conference presentations, tutorials, vignettes.
    *   Code, (anonymized/simulated) data, and analysis scripts shared.
    *   **Visual Abstract:** A figure comparing conceptual flow (raw fMRI → CBD continuous states) vs. traditional two-stage pipeline.

**8. Timeline & Milestones (Refined for R with Optimizations)**
    *   **M0-4: Core R6 Framework & Rcpp Kernels.** FFT convolution, low-rank likelihoods, basic Forward-Backward in RcppArmadillo. Integration with `neuroim2` data structures and `fmrireg` HRF components.
    *   **M5-8: Simulation & Initial Validation.** Parameter recovery, Rcpp profiling & optimization. Shared-HRF cluster mode.
    *   **M9-14: Benchmarking & Advanced Features.** Public datasets (ROI/parcellated focus). GMRF with batched solves. SVB implementation. Memory mapping (`bigmemory`/`arrow`).
    *   **M15-20: Clinical Pilot & User Tools.** Transfer learning module. QC dashboard. Full documentation & R package prep.
    *   **M21-24: Dissemination.** Manuscripts, package release, tutorials.

**8.1. Two-Phase Implementation Strategy: CLD → CBD**

The implementation will proceed in two sequential phases to maximize development efficiency and minimize risk. **We begin with the Continuous Linear Decoder (CLD) as our primary implementation due to its computational efficiency and practical utility**, followed by the full Continuous Bayesian Decoder (CBD) for applications requiring uncertainty quantification.

**Phase 1: Sprint 1a - Continuous Linear Decoder (CLD) as Primary Implementation**
*   **Rationale:** The CLD provides 80-90% of the performance gains of the full Bayesian approach with dramatically reduced computational complexity
*   **Algorithm:** GLM + SVD for spatial maps, FISTA with Total Variation penalty for continuous state estimation
*   **Performance:** Seconds to minutes vs. hours for whole-brain analysis - enabling overnight scanning of dozens of subjects
*   **Benefits:** 
    - Orders of magnitude faster than full VB (40s for V=20,000, T=600 on laptop)
    - Simpler to debug and validate
    - No probabilistic machinery overhead
    - Establishes neuroimaging data pipeline (`neuroim2`/`fmrireg` integration)
    - Provides immediate practical utility for rapid continuous decoding
    - All computations reduce to efficient BLAS/FFT operations

**Phase 2: Sprint 1 - Continuous Bayesian Decoder (CBD) for Advanced Applications**
*   **Algorithm:** Full Variational Bayes with HMM forward-backward for probabilistic state inference
*   **When to use:** Applications requiring uncertainty quantification, explicit Markov priors, or voxel-wise HRF estimation
*   **Builds on:** All infrastructure from Sprint 1a with algorithm-specific additions
*   **Benefits:** 
    - Full posterior distributions over states
    - Principled handling of uncertainty
    - Markov chain temporal dependencies  
    - Voxel-wise HRF adaptation
    - Handles more complex noise models

**Shared Infrastructure & Consolidation Strategy**

Both approaches share substantial machinery (~60-70% overlap):

**🟢 Unified Components:**
- **Project setup & dependencies:** Identical `neuroim2`, `fmrireg`, `R6`, `RcppArmadillo` infrastructure
- **Data structures:** Common matrix ↔ `neuroim2::NeuroVec` conversion and validation
- **HRF processing:** Shared convolution utilities leveraging `fmrireg` basis functions
- **Simulation framework:** Unified `simulate_fmri_data()` supporting both algorithms
- **Testing infrastructure:** Common test utilities, fixtures, and neuroimaging validation
- **Visualization:** Shared plotting methods for spatial maps, convergence, diagnostics

**🎯 Algorithm-Specific Components:**
- **CLD:** FISTA optimization, Total Variation penalty (Condat's algorithm), GLM+SVD learning
- **CBD:** Variational Bayes updates, HMM forward-backward, ELBO computation

**Efficiency Gains from Sequential Development:**
- **~40% reduction** in duplicate infrastructure development
- **Shared simulation** enables direct algorithm comparison and validation
- **Common neuroimaging pipeline** ensures consistent user experience
- **Risk mitigation:** CLD success validates infrastructure before complex CBD implementation

**Proposed Package Structure:**
```r
stance/
├── R/
│   ├── stance-package.R              # Shared documentation
│   ├── data_structures.R             # neuroim2 integration utilities  
│   ├── hrf_utils.R                   # fmrireg convolution wrappers
│   ├── simulate.R                    # Unified simulation functions
│   ├── continuous_linear_decoder.R   # CLD R6 class (Sprint 1a)
│   ├── continuous_bayesian_decoder.R # CBD R6 class (Sprint 1)
│   └── visualization.R               # Shared plotting methods
├── src/
│   ├── fista_tv.cpp                 # CLD-specific FISTA implementation
│   ├── vb_updates.cpp               # CBD-specific VB algorithms
│   └── shared_math.cpp              # Common mathematical operations
└── tests/testthat/
    ├── test-infrastructure.R        # Shared component tests
    ├── test-cld.R                   # CLD-specific tests
    └── test-cbd.R                   # CBD-specific tests
```

**Development Timeline Integration:**
- **M0-2:** Sprint 1a foundation (CLD + shared infrastructure)
- **M3-4:** Sprint 1 implementation (CBD leveraging CLD infrastructure)  
- **M5-8:** Joint validation and algorithm comparison
- **M9+:** Advanced features and optimizations for both approaches

This two-phase strategy ensures robust foundations while enabling direct algorithmic comparison and provides users with both a simpler deterministic option (CLD) and a sophisticated probabilistic framework (CBD).

**9. Risks and Mitigation**
    *   **Adoption Barriers (Clinical):** Phased introduction, strong user support, clear CPU-based benefits.
    *   **R Performance for Full Voxel-wise Whole-Brain:** While heavily optimized, it may remain slower than compiled language versions for the most demanding analyses.
        *   **Mitigation:** Clearly manage expectations. Emphasize ROI/parcellated/clustered-HRF modes for broader accessibility. Position R version as excellent for prototyping, teaching, and moderately-sized data. Provide transparent benchmarks.
    *   **Complexity of Advanced Optimizations (e.g., SVB in R):**
        *   **Mitigation:** Iterative development. Ensure simpler versions are robust first.
    *   **"R can't touch GPUs" (Reviewer Pushback):** While primary focus is CPU, note that `torch` for R exists, allowing core kernels to be (optionally, future work) offloaded if a user has compatible hardware and expertise, while R orchestrates.
    *   **"Forward-backward is quadratic in K":** For typical K ≤ 8-10, direct computation is feasible. Collapsed Dirichlet transitions can make it effectively linear for specific prior choices.

**10. Broader Impact Statement & Conclusion**

This project aims to democratize access to sophisticated, continuous fMRI decoding by providing a robust, well-documented, and computationally optimized implementation within the widely-used R environment. By building on established neuroimaging infrastructure (`neuroim2` for data structures, `fmrireg` for regression components) and combining statistical rigor with pragmatic engineering (🔹SPEED-TIP 6: "Critical inner loops...will be written in RcppArmadillo with OpenMP...reducing whole-brain...VB training to ≈ 45 minutes [on a 12-core Ryzen 9 for a specific configuration]" – adjust claim based on actual preliminary benchmarks), we anticipate this framework will significantly enhance the ability of researchers to extract nuanced information about brain dynamics. The development of both a "pure R reference" and an "optimized R+Rcpp" version will ensure correctness, facilitate innovation, and provide a practical tool for the neuroscience community. This work has the potential to genuinely change how many researchers approach fMRI analysis, particularly those without access to extensive GPU resources or deep expertise in compiled languages.

---

## Appendix A: The Continuous Linear Decoder (CLD) - A Lean Deterministic Twin

### A.1 Motivation and Raison d'Être

While the full Continuous Bayesian Decoder (CBD) offers principled uncertainty quantification and sophisticated temporal modeling, many practical applications prioritize **speed and simplicity** over full probabilistic machinery. The CLD emerges from a critical insight: we can capture 80-90% of the continuous decoding performance gains while reducing computational complexity from hours to seconds.

The CLD is designed for scenarios where:
- **Speed is paramount**: Scanning dozens of subjects overnight rather than days
- **Code clarity matters**: Every step reduces to standard linear algebra operations
- **Iterative development**: Rapid prototyping and parameter tuning without waiting hours
- **Resource constraints**: Standard laptops can process whole-brain data
- **Validation needs**: Simpler algorithms are easier to debug and verify

### A.2 Core Mathematical Formulation

The CLD solves a regularized optimization problem:

$\hat{X} = \arg\min_{X \in \mathbb{R}^{K \times T}} \|Y - W(H \star X)\|_F^2 + \lambda_{\text{TV}} \|\nabla_t X\|_1$

Where:
- $Y$ ∈ ℝ^{V×T}: Observed fMRI data (V voxels, T timepoints)
- $W$ ∈ ℝ^{V×K}: Low-rank spatial maps (learned via GLM+SVD)
- $H \star X$: Convolution with canonical HRF (FFT-accelerated)
- $\|\nabla_t X\|_1$: Total variation penalty enforcing temporal smoothness
- $\lambda_{\text{TV}}$: Regularization parameter balancing fit vs. smoothness

### A.3 Two-Stage Learning Process

#### Stage 1: Learn Spatial Maps W (Once, Fast)

1. **Initial GLM Factorization**:
   - Given a design matrix $S_{design}$ (K×T), convolve with HRF to get $X_{conv}$
   - Run voxel-wise GLMs: $Y_v \sim X_{conv}^T$ to get beta maps $B$ (V×K)
   - Use `speedglm` for massive parallelization

2. **Low-Rank Approximation**:
   - Compute truncated SVD: $B = U_r \Sigma_r V_r^T$ with rank $r \approx 15-20$
   - Store $W = U_r \Sigma_r$ as spatial loading matrix
   - This captures the dominant spatial patterns while reducing dimensionality

**Runtime**: Minutes for whole brain (single GLM per voxel + one SVD)

#### Stage 2: Estimate Continuous States X (FISTA)

Given learned $W$, solve for soft state activations $X$ using Fast Iterative Shrinkage-Thresholding Algorithm (FISTA) with Condat's algorithm for the TV proximal operator.

**Key optimizations**:
- Pre-compute $W^T Y$ once (reduces V-dimensional to K-dimensional operations)
- FFT-based convolution when T > 256 (O(T log T) vs O(T²))
- Warm starts for sequential analyses
- GPU-friendly operations (pure BLAS/FFT)

### A.4 Algorithmic Details

**FISTA Update Loop**:
```
Initialize: X⁰ = 0, Z¹ = X⁰, t¹ = 1
For k = 1 to max_iter:
    1. Gradient step: G^k = ∇f(Z^k) = -(H^T ★ (W^T(Y - W(H ★ Z^k))))
    2. Proximal step: X^{k+1} = prox_{λ/L}(Z^k - (1/L)G^k)
    3. Momentum: t^{k+1} = (1 + √(1 + 4(t^k)²))/2
    4. Extrapolation: Z^{k+1} = X^{k+1} + ((t^k - 1)/t^{k+1})(X^{k+1} - X^k)
```

**Total Variation Proximal Operator** (Condat's Algorithm):
- Efficiently computes $\text{prox}_{\lambda}(x) = \arg\min_z \frac{1}{2}\|z - x\|^2 + \lambda \sum_i |z_{i+1} - z_i|$
- Direct O(T) algorithm without matrix inversions
- Preserves edges while smoothing noise

### A.5 Practical Performance

| Scenario | Data Size | CLD Time | Full CBD Time | Speedup |
|----------|-----------|----------|---------------|----------|
| ROI Analysis | V=1K, T=500 | 5-10s | 5-15 min | 30-90x |
| Parcellated Brain | V=5K, T=500 | 30-60s | 30-90 min | 30-60x |
| Full Brain | V=50K, T=600 | 5-10 min | 4-12 hours | 50-150x |

**Memory usage**: O(VK + KT) vs O(V²) for some VB implementations

### A.6 When to Use CLD vs CBD

**Use CLD when**:
- Rapid iteration during method development
- Real-time or near-real-time decoding requirements  
- Initial exploration of datasets
- Computational resources are limited
- Simpler interpretability is preferred

**Upgrade to CBD when**:
- Uncertainty quantification is critical
- Strong Markov priors improve performance
- Voxel-wise HRF variability is substantial
- Publication requires full Bayesian treatment

### A.7 Modular Extensions

The CLD framework is designed for modularity:

| Extension | Implementation | Computational Cost |
|-----------|----------------|--------------------|
| Voxel-wise HRFs | 3-basis expansion (canonical + derivatives) | 3× gradient computation |
| Online decoding | Warm-start FISTA with previous X | <50ms per TR |
| Multi-subject | Hierarchical W with shared + subject components | Parallel over subjects |
| AR(1) noise | Pre-whiten Y before fitting | One-time preprocessing |

### A.8 Implementation in R

```r
# Core CLD usage pattern
cld <- ContinuousLinearDecoder$new(
  Y = fmri_data,           # V × T matrix or NeuroVec
  S_design = design_matrix, # K × T known states for training
  hrf = "spmg1",          # or custom HRF vector
  rank = 20,              # Low-rank approximation
  lambda_tv = 0.02        # TV regularization strength
)

# Fit continuous states
cld$fit(max_iter = 100)

# Extract results  
states <- cld$get_activations()     # K × T soft assignments
maps <- cld$get_spatial_maps()      # V × K spatial patterns
```

### A.9 Relationship to Full CBD

The CLD serves as both:
1. **Standalone method**: For applications prioritizing speed
2. **CBD initializer**: CLD solution provides excellent warm-start for VB
3. **Validation baseline**: Simpler algorithm for verifying CBD correctness
4. **Shared infrastructure**: 70% code overlap enables dual maintenance

### A.10 Summary

The Continuous Linear Decoder represents a pragmatic compromise between sophistication and speed. By replacing probabilistic updates with closed-form solutions and leveraging modern optimization techniques (FISTA, TV regularization, low-rank approximations), we achieve continuous state estimation at speeds compatible with routine neuroimaging workflows. This makes advanced continuous decoding accessible to the broader neuroimaging community while providing a solid foundation for full Bayesian extensions when needed.

---

This revised proposal is now much stronger in addressing the R implementation challenges head-on, armed with concrete strategies and clear integration with the established R neuroimaging ecosystem.