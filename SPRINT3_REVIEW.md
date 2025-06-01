# Sprint 3 Code Review & Planning

Sprint 3 introduced spatial smoothing of HRF coefficients using a Gaussian Markov Random Field (GMRF) prior and finalised the low‑rank optimisation of likelihood calculations. Optional stochastic variational Bayes (SVB) updates and advanced diagnostics were also added. The package continues to build on `neuroim2` and `fmrireg` for spatial operations and HRF generation.

## Highlights

- **GMRF HRF smoothing** via `create_gmrf_laplacian` and batched solves
- **Low‑rank likelihood computation** fully integrated with Rcpp kernels
- **SVB option** for scalability on larger datasets
- **Diagnostic vignettes** covering spatial priors and performance tuning

## Refactor Notes

1. **Spatial utilities** consolidate neighbour extraction and Laplacian creation.
2. **Profiling helpers** now check proposal performance targets.
3. **Integration tests** cover Sprint 1–3 features with configurable backends.

## Sprint 4 Planning & Risk Assessment

Sprint 4 will apply the model to real fMRI data and polish the user-facing API.
Key tasks include:

- Real‑data ingestion compatible with BIDS/fMRIPrep outputs
- High‑level wrappers and diagnostic dashboards
- Memory‑mapped file support for whole‑brain analysis
- Comprehensive documentation and benchmarking on real datasets

**Risks and Mitigations**

| Risk | Impact | Mitigation |
|------|--------|-----------|
| GMRF solver instability in disconnected regions | High | Detect connected components and solve separately |
| Memory explosion with 50k+ voxels | High | Use out‑of‑core solving and approximate methods |
| SVB convergence issues | Medium | Keep SVB optional, add convergence diagnostics |
| Sparse solver platform failures | Medium | Provide multiple backends and test on all CI platforms |
| Performance regressions | High | Continuous benchmarking with feature flags |

With these plans in place Sprint 4 will focus on validating the full pipeline on real neuroimaging data while keeping runtimes within the proposal targets.
