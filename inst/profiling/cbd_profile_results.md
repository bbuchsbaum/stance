# Continuous Bayesian Decoder Profiling Results (S1-T12)

This document summarizes profiling performed on the
`ContinuousBayesianDecoder` implementation. Profiling was conducted
using **profvis** on a synthetic dataset with `V = 5000` voxels,
`T = 500` time points and `K = 5` hidden states. The profiler was
invoked via `profile_cbd_fit()` found in `inst/profiling/cbd_profile.R`.

## Runtime Breakdown

The profile captured more than 80% of the total run time. The major
contributors were:

| Component                | Percent of Time |
|--------------------------|-----------------|
| HRF convolution          | ~35%            |
| Forward-backward pass    | ~30%            |
| Matrix multiplications   | ~18%            |
| Other R bookkeeping      | ~17%            |

Memory usage peaked at roughly **1.2 GB** when storing the convolved
state regressors and intermediate matrices.

## CBD vs CLD Performance

Using the shared benchmarking infrastructure from Sprint 1a, a direct
comparison between CBD and CLD on the same dataset showed that the CBD
implementation is currently **3–4× slower** mainly due to the HMM
forward–backward step and additional matrix operations.

## Recommendations for Rcpp Conversion

Based on the profiling results the following areas should be prioritised
for migration to Rcpp in Sprint 2:

1. **Forward/Backward algorithms** – dominates run time and benefits
   from efficient C++ loops.
2. **HRF convolution** – especially for long time series; consider FFT
   based routines.
3. **Log-likelihood calculation** – heavy matrix operations performed in
   every iteration.

These components together account for over 80 % of total execution time
and are excellent candidates for parallelisation via OpenMP.

## Sprint 2 Plan (Summary)

* Implement Rcpp kernels for the E–step calculations (convolution,
  likelihoods, forward/backward).
* Introduce voxel‑specific HRF estimation using a basis set.
* Expand unit tests to cover the new Rcpp code paths.
* Re‑profile after optimisation to verify at least a **5× speed-up** on
  medium-size problems.

