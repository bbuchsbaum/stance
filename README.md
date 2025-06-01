# stance

`stance` provides continuous decoding tools for fMRI analysis. Two algorithms are available:

* **Continuous Linear Decoder (CLD)** – deterministic GLM+SVD approach with FISTA optimization.
* **Continuous Bayesian Decoder (CBD)** – probabilistic model using a hidden Markov prior and variational Bayes.

Use **CLD** when you have a reliable design matrix and want fast estimation. Use **CBD** when state boundaries are unknown or you need uncertainty estimates for latent states. See `vignette("getting-started", package = "stance")` for worked examples.

## Parallel Processing

Several low-level routines such as row-wise convolution and total variation denoising are implemented in C++ with optional OpenMP support. When compiled with OpenMP, these functions can run in parallel across rows and the number of threads can be controlled via the `n_threads` argument (use `0` to auto-detect).

## HRF estimation and Rcpp engine

Hemodynamic response functions (HRFs) are generated using helper
functions such as `create_hrf_basis_canonical()` and
`create_hrf_basis_fir()` which internally call the *fmrireg* package.
The resulting basis matrices are orthonormal by default so that unit
tests remain stable.  Voxel‑wise HRF coefficients are initialised using
the Rcpp routine `update_hrf_coefficients_batched_cpp` and stored for
each voxel.

The core E‑step kernels—`compute_log_likelihoods_rcpp`,
`forward_pass_rcpp` and `backward_pass_rcpp`—are implemented with
RcppArmadillo.  When the package is compiled without these functions,
pure R fallbacks are used automatically.  All routines operate on plain
matrices or on `neuroim2::NeuroVec` objects (the data are extracted with
metadata preserved).

## Spatial Priors and Performance Tips

Sprint 3 adds optional GMRF smoothing of voxel-wise HRF coefficients and stochastic VB updates for large datasets. Enable these features when constructing `ContinuousBayesianDecoder` via the `use_gmrf` and `batch_size` arguments. See `vignette("spatial-priors-with-stance", package = "stance")` for worked examples and the `performance-optimization-guide` vignette for detailed SPEED-TIPS and a small decision tree on choosing parameters.

