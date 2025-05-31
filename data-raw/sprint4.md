Okay, Sprint 4 is the final sprint focused on core development before moving towards broader application, packaging, and dissemination. The goals for Sprint 4 are:
1.  **Application to Real fMRI Data:** Test the full model on actual, preprocessed fMRI datasets.
2.  **User-Facing Tools & Polish:** Develop wrapper functions, comprehensive diagnostic dashboards, and improve overall usability.
3.  **Final Performance Optimizations & Robustness:** Address any lingering performance issues, enhance numerical stability, and implement memory management strategies like memory-mapped files if not already done.
4.  **Complete Documentation & Prepare for Release:** Finalize all documentation, vignettes, and prepare the R package for a preliminary release (e.g., on GitHub).
5.  **Benchmarking and Comparative Analysis:** Conduct initial benchmarks against simpler methods using real or highly realistic simulated data.

**Sprint 4 Goal:** Validate the CBD model on real fMRI data. Finalize R package structure, documentation, and user-facing tools. Implement final performance optimizations and robustness checks. Conduct initial benchmarking. Ensure seamless integration with established R neuroimaging infrastructure (`neuroim2`, `fmrireg`) for production-ready deployment. **This sprint builds on all previous work: Sprint 1a shared infrastructure, Sprint 1 CBD implementation, Sprint 2 Rcpp optimizations, and Sprint 3 GMRF/advanced features. Validates that proposal.md performance targets are met in real-world scenarios (ROI: 5-15 min, Parcellated: 30-90 min, Full WB: 4-12 hours).** **All code merged in Sprint 4 must pass `R CMD check --as-cran` (≤ 1 NOTE) on macOS, Linux, and Windows CI before PR approval.**

**Assumptions for Sprint 4:**
*   **All previous sprints completed:** Shared infrastructure (1a), CBD core (1), Rcpp optimizations (2), GMRF and advanced features (3) are functional and tested
*   Access to preprocessed real fMRI datasets (e.g., from HCP, CamCAN, or local data in BIDS format).
*   `neuroim2` and `fmrireg` integration from previous sprints is robust and ready for real-world application.
*   Performance targets from proposal.md have been verified on simulated data
*   **Focus on real data validation**: This sprint transitions from simulation to production use
*   **CI gate:** Matrix build (`macOS-latest`, `ubuntu-latest`, `windows-latest`) must stay green after each ticket merge.

---

**Sprint 4 Tickets:**

**Theme 1: Real fMRI Data Application & Validation**

1.  **Ticket S4-T01: Develop Data Ingestion for Real fMRI Data (BIDS-aware)** 🟢 **[LEVERAGES ALL SPRINTS]**
    *   Description: Create R functions to load preprocessed fMRI data (e.g., NIfTI files) and corresponding event/design files (e.g., BIDS `.tsv` files). This should handle:
        *   Reading 4D NIfTI images using `neuroim2::read_vol()` for seamless integration with `NeuroVol` and `NeuroVec` structures.
        *   **🔄 Add lazy loading:** Use `neuroim2::mmap_vol()` for memory-mapped access to large datasets (🔹SPEED-TIP #8 from proposal)
        *   Extracting voxel time series from specified ROIs or whole-brain masks using `neuroim2` spatial indexing capabilities.
        *   Parsing event files to create initial state indicator sequences for supervised training or to define "rest" periods.
        *   Handling potential issues like different TRs, scan durations, slice timing corrections.
        *   Converting between BIDS/NIfTI formats and `neuroim2` data structures efficiently.
        *   **Add `validate_bids_run()` utility** that checks TR length vs. events and prints warnings for dimension mismatches.
        *   Integration with fMRIPrep outputs (leveraging Sprint 1a data structure utilities)
    *   Acceptance Criteria: Can load preprocessed fMRI data and associated design information into the required R structures for the CBD model, with native `neuroim2` data structure support and robust BIDS validation. **Compatible with all previous sprint components.**
2.  **Ticket S4-T02: Initial Model Fit on Real fMRI Dataset (Single Subject ROI)** 🟡 **[VALIDATES SPRINTS 1-3]**
    *   Description: Apply the full CBD model to a manageable real fMRI dataset (e.g., data from one subject, one task run, focusing on a specific ROI known to be involved in the task).
        *   Carefully select initial parameters (`r`, `K`, `lambda_h`, HRF basis using `fmrireg` HRF components from Sprint 1a).
        *   **🔄 Include motion-regressor option:** Read `_confounds.tsv` files to handle nuisance regressors for day-one demos.
        *   Monitor ELBO convergence (Sprint 1), inspect estimated parameters (`U,V,H_v,Π`), and latent states (`S_gamma`).
        *   Verify Rcpp optimizations from Sprint 2 provide expected speedup
        *   Test GMRF spatial smoothing from Sprint 3 if ROI spans multiple regions
        *   Ensure the model works seamlessly with real `neuroim2::NeuroVec` time series data.
        *   **Performance target from proposal: Single-subject ROI (1000 voxels) completes in 5-15 minutes**
    *   Acceptance Criteria: Model runs to completion on real data. ELBO converges. Estimated parameters and states appear plausible for the given task/ROI. **Meets proposal.md performance target.** All sprint features integrated successfully.
3.  **Ticket S4-T03: Expand to Whole-Brain Analysis (Parcellated)** 🟡 **[STRESS-TESTS ALL OPTIMIZATIONS]**
    *   Description: Scale up to parcellated whole-brain analysis per proposal strategy:
        *   **Primary approach - Parcellation:** Use a brain atlas (e.g., AAL, Schaefer) to define ~300-400 parcels. Average BOLD signal within parcels using `neuroim2` ROI functionality. Run CBD on parcel time series (~5000 effective voxels).
        *   Leverage all optimizations: low-rank projection (Sprint 2/3), FFT convolution (Sprint 2/3), GMRF (Sprint 3)
        *   Option 2 (Clustered HRF): Use the "clustered-HRF" mode (from Sprint 3) where HRFs are estimated per spatial cluster
        *   **Performance target from proposal: Parcellated whole-brain completes in 30-90 minutes**
        *   **Memory target from proposal: <8GB RAM usage**
        *   Document full voxel-wise requirements: ≥64 GB RAM, 4-12 hours runtime
    *   Acceptance Criteria: Model runs on parcellated whole-brain data **within proposal.md timeframe**. Results are interpretable. Memory usage stays within bounds. Sprint 3 optimizations demonstrably effective.
4.  **Ticket S4-T04: Qualitative Validation of Results on Real Data** 🟢 **[USES ALL INFRASTRUCTURE]**
    *   Description: For the real data analyses, leverage all sprint components:
        *   Compare estimated spatial maps `W_k` with known functional neuroanatomy using Sprint 3 visualization tools
        *   Examine estimated HRF shapes using Sprint 1a `fmrireg` HRF visualization
        *   Correlate inferred state probability time courses `S_gamma` with experimental design
        *   Validate GMRF smoothing (Sprint 3) produces anatomically plausible spatial patterns
        *   Use Sprint 3 diagnostic tools to assess model fit quality
        *   **Auto-open in neuroimaging viewers:** Export results as NIfTI using `neuroim2::write_vol()`
        *   Compare performance metrics against baseline two-stage approach
    *   Acceptance Criteria: Model outputs show meaningful correspondence with neuroscience knowledge. **All visualization and diagnostic tools from previous sprints function correctly with real data.** Clear demonstration of CBD advantages.

**Theme 2: User Interface, Diagnostics & Package Polish**

5.  **Ticket S4-T05: Develop High-Level Wrapper Functions** 🟡 **[INTEGRATES ALL COMPONENTS]**
    *   Description: Create user-friendly wrapper functions that abstract complexity:
        ```r
        # Main wrapper function
        run_cbd_analysis <- function(
          data,                    # Matrix or NeuroVec/NeuroVol
          n_states = 3,           
          rank = 10,
          algorithm = c("CBD", "CLD"),  # Support both from Sprint 1a/1
          use_gmrf = TRUE,        # Sprint 3 feature
          lambda_h = 1.0,
          hrf_basis = fmrireg::HRF_SPMG1,  # Sprint 1a integration
          optimization = c("full", "fast", "memory_efficient"),
          verbose = TRUE,
          ...
        ) {
          # Auto-detect data type using Sprint 1a utilities
          # Select appropriate optimizations based on data size
          # Configure based on optimization profile
          # Run appropriate algorithm
          # Return standardized results object
        }
        ```
        *   Support both CBD (Sprint 1) and CLD (Sprint 1a) algorithms
        *   Auto-configure based on data size and available resources
        *   Provide sensible defaults based on Sprint 1-4 experience
        *   Return unified results object with S3 methods
    *   Acceptance Criteria: Users can run analyses with minimal code. **Seamlessly integrates all sprint features.** Clear parameter documentation. Intelligent defaults.
6.  **Ticket S4-T06: Comprehensive Diagnostic Dashboard** 🟢 **[EXTENDS SPRINT 3 DIAGNOSTICS]**
    *   Description: Create production-ready diagnostic reporting building on Sprint 3 work:
        *   Extend Sprint 3 diagnostic functions to generate full HTML reports
        *   Include all diagnostic plots from Sprints 1-3 in organized sections
        *   Add new real-data specific diagnostics:
            - Motion parameter correlations
            - Comparison with anatomical atlases
            - Cross-validation results if available
        *   Performance metrics section showing:
            - Runtime vs. proposal targets
            - Memory usage vs. limits
            - Optimization effectiveness (with/without features)
        *   **Interactive elements:** Use `plotly` for key plots
        *   **Export functionality:** Save key results as .RData, CSV, NIfTI
        *   **Target: Report generation <2 min, HTML size <10MB**
    *   Acceptance Criteria: Comprehensive diagnostic report integrating all sprint features. **Performance metrics clearly shown against proposal targets.** Publication-ready figures. Fast generation.
7.  **Ticket S4-T07: Finalize Package Structure** 🟡 **[FINAL INTEGRATION]**
    *   Description: Complete R package structure incorporating all sprints:
        ```
        stance/
        ├── DESCRIPTION         # All dependencies from Sprints 1a-4
        ├── NAMESPACE          # Proper exports/imports
        ├── NEWS.md            # Chronicle all sprint deliverables
        ├── R/
        │   ├── stance-package.R (Sprint 1a)
        │   ├── data_structures.R (Sprint 1a)
        │   ├── hrf_utils.R (Sprint 1a)
        │   ├── simulate.R (Sprint 1a/1)
        │   ├── continuous_linear_decoder.R (Sprint 1a)
        │   ├── continuous_bayesian_decoder.R (Sprint 1)
        │   ├── visualization.R (Sprint 1a/3)
        │   ├── spatial_utils.R (Sprint 3)
        │   ├── cross_validation.R (Sprint 3)
        │   ├── wrapper_functions.R (Sprint 4)
        │   └── diagnostics.R (Sprint 3/4)
        ├── src/
        │   ├── fista_tv.cpp (Sprint 1a)
        │   ├── vb_updates.cpp (Sprint 1/2)
        │   ├── likelihood_lowrank.cpp (Sprint 2/3)
        │   ├── convolution_fft.cpp (Sprint 2/3)
        │   └── gmrf_batch_solver.cpp (Sprint 3)
        └── vignettes/
            ├── getting-started.Rmd
            ├── performance-guide.Rmd
            └── real-data-analysis.Rmd
        ```
    *   Ensure `R CMD check --as-cran` passes with ≤1 NOTE
    *   Configure `pkgdown` site with all documentation
    *   Acceptance Criteria: **Complete integration of all sprint components.** Package structure follows best practices. Ready for release.
8.  **Ticket S4-T08: Comprehensive Vignettes** 🟢 **[DOCUMENTS ALL FEATURES]**
    *   Description: Create vignettes showcasing the complete system:
        *   **"Getting Started"**: Simple example using Sprint 1a simulation, basic CBD (Sprint 1)
        *   **"Performance Optimization"**: Demonstrate Sprint 2/3 optimizations, when to use each
        *   **"Real Data Analysis"**: Complete workflow from BIDS to results using all features
        *   **"Algorithm Comparison"**: CLD vs CBD, when to use each
        *   Include performance benchmarks showing proposal target achievement
        *   Use consistent neuroimaging data throughout
        *   **Include mini BIDS dataset** (~50MB) for reproducibility
    *   Acceptance Criteria: Vignettes demonstrate all sprint features. **Clear evidence of meeting proposal performance targets.** Reproducible examples.
9.  **Ticket S4-T09: Complete Documentation** 🟡 **[FINAL POLISH]**
    *   Description: Ensure comprehensive documentation:
        *   All functions from Sprints 1a-4 have complete roxygen2 docs
        *   Examples show both matrix and neuroimaging data usage
        *   Performance tips based on Sprint 2/3 optimizations documented
        *   Troubleshooting guide for common issues discovered during testing
        *   **API stability commitment**: Mark stable vs experimental features
        *   Cross-references between related functions
    *   Acceptance Criteria: No undocumented objects. Examples run successfully. **Clear guidance on using optimizations to meet performance targets.**

**Theme 3: Performance Validation & Robustness**

10. **Ticket S4-T10: Memory-Mapped File Handling** 🟡 **[COMPLETES 🔹SPEED-TIP #8]**
    *   Description: Finalize out-of-memory support for large datasets:
        *   Complete integration with `neuroim2` memory-mapped volumes
        *   Implement `arrow`-based solution for time series data
        *   Windows fallback using `bigmemory::filebacked.big.matrix`
        *   Automatic switching based on available RAM vs data size
        *   Benchmark memory-mapped vs in-memory performance
        *   **Target: Handle 100GB+ datasets on 16GB machines**
    *   Acceptance Criteria: Large datasets processable without memory errors. **Enables full whole-brain analysis per proposal.** Cross-platform compatibility.
11. **Ticket S4-T11: Final Performance Validation** 🔴 **[CRITICAL - VALIDATES PROPOSAL]**
    *   Description: Comprehensive performance validation against proposal targets:
        ```r
        validate_performance_targets <- function() {
          scenarios <- list(
            roi = list(
              voxels = 1000, T = 500, K = 3, r = 10,
              target_time = c(5, 15),  # minutes
              target_ram = 2  # GB
            ),
            parcellated = list(
              voxels = 5000, T = 500, K = 4, r = 15,
              target_time = c(30, 90),  # minutes
              target_ram = 8  # GB
            ),
            full_wb = list(
              voxels = 50000, T = 1000, K = 5, r = 20,
              target_time = c(240, 720),  # minutes (4-12 hours)
              target_ram = 32  # GB
            )
          )
          
          # Test each scenario with all optimizations
          # Record actual vs target performance
          # Generate performance report
        }
        ```
    *   Document any gaps between targets and actual performance
    *   Provide optimization recommendations for each scenario
    *   Acceptance Criteria: **Clear pass/fail on each proposal.md performance target.** Documented optimization strategies for meeting targets.
12. **Ticket S4-T12: Robustness & Error Handling** 🟡 **[PRODUCTION HARDENING]**
    *   Description: Ensure production-ready robustness:
        *   Add comprehensive input validation using Sprint 1a utilities
        *   Handle edge cases discovered during real data testing:
            - Missing timepoints
            - Motion outliers
            - Singular covariance matrices
            - Disconnected brain regions
        *   Informative error messages referencing neuroimaging concepts
        *   Recovery strategies for numerical issues
        *   Logging framework for debugging
        *   **Add stress tests**: Pathological inputs, extreme parameters
    *   Acceptance Criteria: Graceful handling of all common issues. **Production-ready error handling.** Clear diagnostic messages.

**Theme 4: Benchmarking & Comparative Analysis**

13. **Ticket S4-T13: Standardized Benchmarking Framework** 🟢 **[DEMONSTRATES VALUE]**
    *   Description: Create comprehensive benchmarking comparing CBD to standard approaches:
        *   Implement baseline: `fmrireg::glm_hrf()` + classification (two-stage MVPA)
        *   Use consistent data handling via `neuroim2` for fair comparison
        *   Metrics to collect:
            - Decoding accuracy (AUC, F1, balanced accuracy)
            - Temporal precision (state transition detection)
            - Computational resources (time, memory)
            - HRF recovery accuracy (if ground truth available)
        *   Support both simulated and real data benchmarks
        *   Generate publication-ready comparison tables and figures
    *   Acceptance Criteria: Fair comparison framework. **CBD advantages clearly demonstrated.** Results suitable for manuscript.
14. **Ticket S4-T14: Real Data Benchmarks** 🟡 **[VALIDATES ON REAL DATA]**
    *   Description: Execute benchmarks on multiple real datasets:
        *   **Dataset 1**: Simple motor task (known ground truth)
        *   **Dataset 2**: Complex cognitive task (multiple states)
        *   **Dataset 3**: Resting state (validate "rest" state detection)
        *   Compare CBD with all features vs baseline approaches
        *   Statistical analysis with bootstrap confidence intervals
        *   Assess impact of each optimization (ablation study)
        *   **Document minimum performance gains needed to justify complexity**
    *   Acceptance Criteria: Benchmarks on ≥3 real datasets. **Statistical evidence of CBD advantages.** Clear recommendations for use cases.
15. **Ticket S4-T15: (Optional) Transfer Learning Demo** 🟢 **[STRETCH GOAL]**
    *   Description: If time permits, demonstrate transfer learning capabilities:
        *   Train on one dataset, fine-tune on another
        *   Show improved convergence and performance
        *   Handle different scanners/protocols via `neuroim2` coordinate systems
        *   Create example showing clinical application potential
    *   Acceptance Criteria: Basic transfer learning demonstrated. **Future research directions identified.**

**Theme 5: Release Preparation & Future Planning**

16. **Ticket S4-T16: Code Quality & Release Prep** 🟡 **[FINAL QUALITY CHECK]**
    *   Description: Ensure release-ready code quality:
        *   Run all linters: `lintr::lint_package()`, `goodpractice::gp()`
        *   Apply consistent style: `styler::style_pkg()`
        *   Security audit: Check for any hardcoded paths, credentials
        *   License compliance: Verify all dependencies compatible
        *   **Performance regression tests**: Ensure no sprint degraded performance
        *   Create release checklist for future versions
    *   Acceptance Criteria: Clean code quality reports. **No performance regressions across sprints.** Release checklist created.
17. **Ticket S4-T17: Project Summary & Demo** 🟢 **[SHOWCASES ACHIEVEMENT]**
    *   Description: Create compelling project summary:
        *   Executive summary of achievements across all sprints
        *   Performance achievements vs proposal targets (table/figures)
        *   Real data analysis showcase with visualizations
        *   **Live demo video** showing complete workflow
        *   Slide deck for conference presentations
        *   **Emphasize integration with R neuroimaging ecosystem**
    *   Acceptance Criteria: Professional project summary. **Clear demonstration of meeting proposal goals.** Ready for dissemination.
18. **Ticket S4-T18: Future Roadmap** 🟢 **[PLANS BEYOND SPRINTS]**
    *   Description: Document future development plans:
        *   **Immediate**: GitHub release, Neuroconductor submission
        *   **Short-term** (3-6 months):
            - Community feedback incorporation
            - Additional HRF basis functions
            - Multi-session/subject extensions
        *   **Medium-term** (6-12 months):
            - GPU acceleration via `torch` for R
            - Advanced HMM variants
            - Clinical validation studies
        *   **Long-term** (1+ years):
            - Integration with other neuroimaging packages
            - Real-time analysis capabilities
            - Clinical deployment guidelines
    *   Acceptance Criteria: Clear development roadmap. **Realistic timeline based on sprint velocity.** Community engagement plan.

---

**Cross-Cutting Risks & Mitigations (Sprint 4):**

| Risk | Impact | Mitigation |
|------|--------|------------|
| Real data doesn't meet assumptions | High | Start with well-validated datasets, provide diagnostic tools |
| Performance targets not met on real data | High | Profile continuously, document optimization strategies |
| Memory requirements exceed targets | High | Default to parcellation, document full-brain requirements |
| Package too complex for users | Medium | Provide simple wrapper functions with smart defaults |
| CI timeouts with real data tests | Medium | Use smaller test datasets, cache preprocessed data |
| Cross-platform issues emerge late | Medium | Test on all platforms throughout sprint |
| Documentation becomes outdated | Low | Automated checks, vignette testing in CI |
| Benchmark results disappointing | Medium | Focus on specific use cases where CBD excels |

---

**Sprint 4 Exit Criteria:**
1. **Functionality**: All features working on real data
   - Real fMRI data pipeline complete and tested
   - All optimizations functioning correctly
   - Package passes CRAN checks (≤1 NOTE)
2. **Performance**: Meeting proposal.md targets on real data
   - ROI: 5-15 minutes ✓
   - Parcellated: 30-90 minutes ✓  
   - Full WB: 4-12 hours (documented) ✓
   - Memory usage within specified bounds ✓
3. **Quality**: Production-ready package
   - Test coverage ≥90%
   - Documentation complete and accurate
   - Vignettes reproducible
   - No performance regressions
4. **Validation**: Demonstrated advantages
   - Benchmarks show CBD benefits
   - Real data results scientifically valid
   - User feedback incorporated
5. **Release Ready**: Package prepared for distribution
   - GitHub release tagged
   - Neuroconductor submission prepared
   - Community engagement plan in place
   - Future roadmap documented

**Sprint 4 Summary:** This final sprint transforms the CBD implementation from a research prototype into a production-ready R package for the neuroimaging community. By validating on real data, finalizing all optimizations, and ensuring seamless integration of all components from Sprints 1a through 3, we deliver a powerful tool that achieves the ambitious performance targets set in the proposal while maintaining scientific rigor and usability. The package is positioned as a significant advancement in fMRI analysis methodology, ready for broad adoption and future development. 🎉