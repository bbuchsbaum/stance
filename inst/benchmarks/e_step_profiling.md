# E-Step Profiling Results (S2-T04)

This report summarises profiling of the **E-step** routines from the
`ContinuousBayesianDecoder`. Profiling targeted the following private
methods:

* `private$.compute_log_likelihoods`
* `private$.forward_backward` (combined forward and backward passes)

Profiling was executed with `profvis`, `lineprof`, and `Rprofmem` on a
synthetic dataset with `V = 2000` voxels, `T = 300` time points and
`K = 4` hidden states. The script used to run the analysis is located in
`inst/profiling/e_step_profile.R`.

## Major Hotspots

| Component                      | Approx. Percent of Time |
|--------------------------------|-------------------------|
| Likelihood state loop          | ~40%                    |
| Forward pass                   | ~30%                    |
| Backward pass                  | ~20%                    |
| Memory allocation (Rprofmem)   | ~10%                    |

The line profiler confirmed that the inner loops within
`compute_log_likelihoods()` and the forward/backward time loops dominate
runtime. Memory profiling revealed transient allocations totalling around
80–100 MB during the forward–backward updates.

## Recommendations

1. Convert the likelihood and forward/backward loops to C++ with
   vectorised operations and possible OpenMP parallelism.
2. Pre-compute reusable terms outside the loops to reduce allocations.
3. Re-run this profiling after the Rcpp implementation to verify speed
   improvements of at least 5×.
